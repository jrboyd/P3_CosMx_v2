"""End-to-end test of reviewer.html in headless Firefox, driven over WebDriver BiDi.

    sbatch --mem=8G --cpus-per-task=2 --time=0:20:00 scripts/sbatch/run_python.sh \
        pathology_review_app/build/test_app_headless.py test

Runs against the BUILT bundle in output/<size>/send/ (reviewer.html + the zip), through the
real file pickers (input.setFiles), and checks what a reader would depend on:

  open zip via the picker -> images decode -> composite/grid render -> keyboard answers ->
  validation blocks an incomplete save -> gate=No skips the form -> revisions append ->
  Export downloads a CSV -> answers survive a page reload -> after storage is wiped,
  Import restores them.

Screenshots land in output/<size>/apptest/ for a human look. Standard library only: no
selenium/playwright on this cluster, so the WebSocket client is the minimal one below.
Firefox is /usr/bin/firefox (ESR 140); BiDi is built in.
"""
import base64
import glob
import json
import os
import shutil
import socket
import struct
import subprocess
import sys
import time

SIZE = sys.argv[1] if len(sys.argv) > 1 else "test"
ROOT = os.path.abspath(os.path.join("pathology_review_app", "output", SIZE))
SEND = os.path.join(ROOT, "send")
T = os.path.join(ROOT, "apptest")
PORT = 20000 + os.getpid() % 20000


class WS:
    def __init__(self, host, port, path):
        self.s = socket.create_connection((host, port), timeout=120)
        key = base64.b64encode(os.urandom(16)).decode()
        self.s.sendall(("GET %s HTTP/1.1\r\nHost: %s:%d\r\nUpgrade: websocket\r\nConnection: Upgrade\r\n"
                        "Sec-WebSocket-Key: %s\r\nSec-WebSocket-Version: 13\r\n\r\n" % (path, host, port, key)).encode())
        buf = b""
        while b"\r\n\r\n" not in buf:
            buf += self.s.recv(4096)
        head, rest = buf.split(b"\r\n\r\n", 1)
        if b" 101 " not in head.split(b"\r\n")[0]:
            raise RuntimeError("websocket upgrade failed: %r" % head[:200])
        self.buf = bytearray(rest)

    def _n(self, n):
        while len(self.buf) < n:
            c = self.s.recv(1 << 20)
            if not c:
                raise EOFError("socket closed")
            self.buf += c
        out = bytes(self.buf[:n])
        del self.buf[:n]
        return out

    def _frame(self, op, payload):
        h = bytearray([0x80 | op])
        n = len(payload)
        if n < 126:
            h.append(0x80 | n)
        elif n < 65536:
            h.append(0x80 | 126); h += struct.pack(">H", n)
        else:
            h.append(0x80 | 127); h += struct.pack(">Q", n)
        m = os.urandom(4)
        h += m
        mm = (m * (n // 4 + 1))[:n]
        self.s.sendall(bytes(h) + (int.from_bytes(payload, "big") ^ int.from_bytes(mm, "big")).to_bytes(n, "big")
                       if n else bytes(h))

    def send(self, text):
        self._frame(0x1, text.encode())

    def recv(self):
        msg = b""
        while True:
            b1, b2 = self._n(2)
            op, n = b1 & 0x0F, b2 & 0x7F
            if n == 126:
                n = struct.unpack(">H", self._n(2))[0]
            elif n == 127:
                n = struct.unpack(">Q", self._n(8))[0]
            mask = self._n(4) if b2 & 0x80 else None
            data = self._n(n)
            if mask:
                data = bytes(b ^ mask[i % 4] for i, b in enumerate(data))
            if op == 0x8:
                raise EOFError("closed by browser")
            if op == 0x9:
                self._frame(0xA, data)
                continue
            msg += data
            if b1 & 0x80:
                return msg.decode()


class Bidi:
    def __init__(self, ws):
        self.ws, self.i = ws, 0

    def cmd(self, method, params):
        self.i += 1
        self.ws.send(json.dumps({"id": self.i, "method": method, "params": params}))
        while True:
            m = json.loads(self.ws.recv())
            if m.get("id") == self.i:
                if m.get("type") == "error":
                    raise RuntimeError("%s: %s %s" % (method, m.get("error"), m.get("message")))
                return m["result"]


def main():
    zips = glob.glob(os.path.join(SEND, "*.zip"))
    assert len(zips) == 1, "expected one bundle zip in %s" % SEND
    zpath = zips[0]
    if os.path.isdir(T):
        shutil.rmtree(T)
    prof, dl = os.path.join(T, "profile"), os.path.join(T, "downloads")
    os.makedirs(prof); os.makedirs(dl)
    with open(os.path.join(prof, "user.js"), "w") as fh:
        for k, v in {"browser.download.folderList": 2, "browser.download.dir": dl,
                     "browser.download.useDownloadDir": True,
                     "browser.download.always_ask_before_handling_new_types": False,
                     "browser.shell.checkDefaultBrowser": False,
                     "datareporting.policy.dataSubmissionEnabled": False}.items():
            fh.write("user_pref(%s, %s);\n" % (json.dumps(k), json.dumps(v)))
    ff = subprocess.Popen(["firefox", "--headless", "--no-remote", "--profile", prof,
                           "--remote-debugging-port", str(PORT)],
                          stdout=open(os.path.join(T, "firefox.log"), "w"), stderr=subprocess.STDOUT)
    results, ok = [], True

    def check(name, cond, detail=""):
        nonlocal ok
        ok &= bool(cond)
        line = "%s  %s%s" % ("PASS" if cond else "FAIL", name, ("  -- " + str(detail)) if detail else "")
        results.append(line)
        print(line, flush=True)

    try:
        for _ in range(120):
            try:
                ws = WS("127.0.0.1", PORT, "/session")
                break
            except (ConnectionRefusedError, OSError):
                time.sleep(0.5)
        else:
            raise RuntimeError("firefox never opened the BiDi port")
        b = Bidi(ws)
        b.cmd("session.new", {"capabilities": {"alwaysMatch": {"unhandledPromptBehavior": {"default": "accept"}}}})
        ctx = b.cmd("browsingContext.getTree", {})["contexts"][0]["context"]
        b.cmd("browsingContext.setViewport", {"context": ctx, "viewport": {"width": 1600, "height": 1000}})
        url = "file://" + os.path.join(SEND, "reviewer.html")

        def ev(expr):
            r = b.cmd("script.evaluate", {"expression": expr, "target": {"context": ctx},
                                         "awaitPromise": True, "resultOwnership": "root"})
            if r["type"] == "exception":
                raise RuntimeError("JS: %s" % json.dumps(r["exceptionDetails"])[:800])
            return r["result"]

        def js(expr):          # expression must evaluate to a JSON-able value
            return json.loads(ev("(async () => JSON.stringify(await (async () => { %s })()))()" % expr)["value"])

        def shot(name):
            d = b.cmd("browsingContext.captureScreenshot", {"context": ctx})["data"]
            with open(os.path.join(T, name), "wb") as fh:
                fh.write(base64.b64decode(d))

        def set_file(elem_id, path):
            node = ev("document.getElementById(%s)" % json.dumps(elem_id))
            b.cmd("input.setFiles", {"context": ctx, "element": {"sharedId": node["sharedId"]}, "files": [path]})

        def wait(cond_js, secs=60):
            for _ in range(secs * 5):
                if js("return !!(%s);" % cond_js):
                    return True
                time.sleep(0.2)
            return False

        def open_bundle():
            b.cmd("browsingContext.navigate", {"context": ctx, "url": url, "wait": "complete"})
            js('document.getElementById("readerName").value = "Test Reader"; return 1;')
            set_file("fileZip", zpath)
            return wait("S.disp && S.disp.comp", 90)

        def key(*ks):
            js("for (const k of %s) document.dispatchEvent(new KeyboardEvent('keydown', {key: k, bubbles: true})); return 1;"
               % json.dumps(list(ks)))

        # ---- 1. start screen, open the zip -------------------------------------------
        b.cmd("browsingContext.navigate", {"context": ctx, "url": url, "wait": "complete"})
        shot("01_start.png")
        js("localStorage.clear(); return 1;")
        check("bundle opens via file picker and first cell renders", open_bundle())
        info = js("return {n: S.manifest.items.length, w: S.disp.data.w, idx: S.idx, "
                  "storage: S.storageOK, err: document.getElementById('banner').classList.contains('hidden') ? '' : "
                  "document.getElementById('banner').textContent};")
        check("manifest items", info["n"] > 0, info)
        check("images decode at window size", info["w"] == 1247, info["w"])
        check("browser storage available on file://", info["storage"], info["err"])
        st = js("const c = S.disp.comp, d = c.getContext('2d').getImageData(0, 0, c.width, c.height).data;"
                "let s = [0, 0, 0]; for (let i = 0; i < d.length; i += 4 * 97) { s[0] += d[i]; s[1] += d[i+1]; s[2] += d[i+2]; }"
                "const n = d.length / (4 * 97); return s.map(v => +(v / n).toFixed(1));")
        check("composite is not blank", sum(st) > 30, "mean RGB %s" % st)
        shot("02_composite.png")
        js("setMode('grid'); return 1;"); shot("03_grid.png")
        js("setMode('DNA'); S.view.z = 4; draw(); return 1;"); shot("04_dna_zoom4.png")
        js("setMode('composite'); fitView(); draw(); return 1;")

        # ---- 2. answer by keyboard ------------------------------------------------------
        nq = js("return S.manifest.form.checklist.length;")
        key("Enter")
        msg = js("return document.getElementById('formMsg').textContent;")
        check("empty form cannot be saved", "Still needed" in msg, msg[:80])
        # gate by Y, then the checklist by digits 1/2/3 (= yes/no/unsure) -- the highlight must
        # walk on by itself and land on the score row after the last question
        key("y", *(["1", "2", "3"] * 3)[:nq])
        a = js("return {act: S.active, conf: ROW.conf(), checks: S.form.checks};")
        exp = (["yes", "no", "unsure"] * 3)[:nq]
        check("digits 1/2/3 answer questions and the highlight walks to the score",
              list(a["checks"].values()) == exp and a["act"] == a["conf"], a)
        key("2")                                   # on the score row: 2 is the score, not "No"
        a = js("return {act: S.active, flag: ROW.flag(), conf: S.form.hrs_confidence};")
        check("on the score row a digit sets the score, highlight moves to flag", a["conf"] == 2 and a["act"] == a["flag"], a)
        key("1")                                   # on the flag row: 1 = flag yes
        key("ArrowDown")
        a = js("return {act: S.active, notes: ROW.notes(), focus: document.activeElement.id, flag: S.form.flag};")
        check("flag row takes 1 = yes; ArrowDown continues into the notes box", a["flag"] is True and
              a["act"] == a["notes"] and a["focus"] == "notes", a)
        js("document.getElementById('notes').blur(); return 1;")
        key("ArrowUp", "ArrowUp", "ArrowUp")
        a = js("return S.active;")
        check("ArrowUp walks back up through score into the questions", a == nq, a)
        key("Escape")
        js("const t = document.getElementById('notes'); t.value = 'note, with \"quotes\"\\nand a newline'; "
           "t.dispatchEvent(new Event('input')); return 1;")
        shot("05_filled.png")
        key("Enter")
        r1 = js("return {idx: S.idx, recs: S.store.records.length, last: S.store.records.at(-1)};")
        check("complete form saves and advances", r1["recs"] == 1 and r1["idx"] == 1, {k: r1[k] for k in ("idx", "recs")})
        check("record carries checklist + score + flag, no identity column", r1["last"]["adequate"] == "yes" and
              r1["last"]["single_cell"] == "yes" and r1["last"]["large_cell"] == "no" and
              r1["last"]["hrs_confidence"] == 2 and r1["last"]["flag"] == "yes" and "identity" not in r1["last"],
              {k: r1["last"].get(k) for k in ("adequate", "single_cell", "large_cell", "hrs_confidence", "flag")})
        res = js("return {resume: !document.getElementById('resumeRow').classList.contains('hidden'), "
                 "fsa: !!window.showOpenFilePicker};")
        check("no resume button where the picker API is absent (Firefox)", not res["fsa"] and not res["resume"], res)
        # item 2: gate No skips the rest
        key("n", "Enter")
        r2 = js("return {idx: S.idx, last: S.store.records.at(-1)};")
        check("gate=No saves with the form skipped", r2["idx"] == 2 and r2["last"]["adequate"] == "no"
              and r2["last"]["hrs_confidence"] == "", r2["last"].get("adequate"))
        # item 3: start, leave unsaved -> draft
        key("y", "y", "ArrowRight")
        d3 = js("return Object.keys(S.store.drafts).length;")
        check("partial answers kept as a draft when skipping", d3 == 1, d3)
        # back to item 1, change the score -> revision 2
        js("goTo(0); return 1;"); wait("S.disp && S.disp.it === S.manifest.items[0]")
        key("4")                                   # score row NOT highlighted: must be ignored
        a = js("return S.form.hrs_confidence;")
        check("a digit is ignored unless the score row is highlighted", a == 2, a)
        js("document.getElementById('rowConf').dispatchEvent(new MouseEvent('mousedown', {bubbles: true})); return 1;")
        key("4", "Enter")
        rv = js("return S.store.records.filter(r => r.item_id === S.manifest.items[0].id).map(r => [r.revision, r.hrs_confidence]);")
        check("editing a saved cell appends revision 2", rv == [[1, 2], [2, 4]], rv)
        shot("06_after_saves.png")

        # ---- 3. export ------------------------------------------------------------------
        csv_text = js("return buildCSV();")
        lines = [l for l in csv_text.lstrip("﻿").split("\r\n") if l]
        # rows end in \r\n; the note's own \n stays inside its quoted field
        check("CSV has header + 3 records (quoted newline inside a field)", len(lines) == 4 and
              'note, with ""quotes""' in csv_text, "%d physical lines" % len(lines))
        js("exportCSV(); return 1;")
        got = []
        for _ in range(50):
            got = [f for f in os.listdir(dl) if f.endswith(".csv")]
            if got:
                break
            time.sleep(0.2)
        check("Export downloads a timestamped CSV", len(got) == 1 and got[0].startswith("HRS_review_Test_Reader_P3HRS-"), got)
        exported = os.path.join(dl, got[0]) if got else None

        # ---- 4. persistence across reload ----------------------------------------------
        ok_open = open_bundle()
        p = js("return {recs: S.store.records.length, drafts: Object.keys(S.store.drafts).length, idx: S.idx};")
        check("answers survive closing and reopening", ok_open and p["recs"] == 3 and p["drafts"] == 1, p)
        check("reopens at the first unsaved cell", p["idx"] == 2, p["idx"])

        # ---- 5. storage wiped -> import restores ---------------------------------------
        js("localStorage.clear(); return 1;")
        open_bundle()
        z = js("return S.store.records.length;")
        set_file("fileImport", exported)
        wait("S.store.records.length === 3", 20)
        q = js("return {recs: S.store.records.length, msg: document.getElementById('formMsg').textContent, "
               "notes: S.store.records.find(r => r.revision == 1 && r.adequate === 'yes').notes};")
        check("import restores answers after storage is wiped", z == 0 and q["recs"] == 3, q["msg"])
        check("notes round-trip intact", q["notes"] == 'note, with "quotes"\nand a newline', repr(q["notes"]))
        shot("07_after_import.png")
        js("document.getElementById('help').classList.remove('hidden'); return 1;"); shot("08_help.png")
    except Exception as e:
        check("test harness ran to completion", False, repr(e))
    finally:
        try:
            ff.terminate(); ff.wait(10)
        except Exception:
            ff.kill()
        shutil.rmtree(prof, ignore_errors=True)
    with open(os.path.join(T, "results.txt"), "w") as fh:
        fh.write("\n".join(results) + "\n%s\n" % ("ALL PASS" if ok else "FAILURES"))
    print("ALL PASS" if ok else "FAILURES")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
