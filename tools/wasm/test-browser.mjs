// npm install --prefix tools/wasm, then: node tools/wasm/test-browser.mjs target/smlib-web-validation
import {createServer} from 'node:http';
import {readFile} from 'node:fs/promises';
import {createHash} from 'node:crypto';
import {resolve, extname, sep} from 'node:path';
const {chromium, firefox, webkit} = await import(process.env.PLAYWRIGHT_MODULE || 'playwright');
const root = resolve(process.argv[2] || 'target/smlib-web-validation');
const manifest = JSON.parse(await readFile(resolve(root, 'manifest.json'), 'utf8'));
for (const [file, expected] of Object.entries(manifest.files)) {
    const actual = createHash('sha256').update(await readFile(resolve(root, file))).digest('hex');
    if (actual !== expected) throw new Error(`Package checksum mismatch: ${file}`);
}
const server = createServer(async (req, res) => {
    try {
        if (req.url === '/') { res.setHeader('Content-Type', 'text/html'); res.end('<!doctype html><title>Wasm production test</title>'); return; }
        const file = resolve(root, '.' + new URL(req.url, 'http://localhost').pathname);
        if (!file.startsWith(root + sep)) { res.writeHead(403).end(); return; }
        const types = {'.wasm':'application/wasm', '.js':'text/javascript', '.mjs':'text/javascript', '.html':'text/html'};
        res.setHeader('Content-Type', types[extname(file)] || 'application/octet-stream');
        res.end(await readFile(file));
    } catch { res.writeHead(404).end(); }
});
await new Promise(done => server.listen(0, '127.0.0.1', done));
let browser;
try {
    const engine = {chromium, firefox, webkit}[process.env.BROWSER || 'chromium'];
    if (!engine) throw new Error('BROWSER must be chromium, firefox or webkit');
    browser = await engine.launch({headless: true, ...(process.env.BROWSER_EXECUTABLE ? {executablePath:process.env.BROWSER_EXECUTABLE} : {})});
    const page = await browser.newPage();
    page.on('pageerror', error => console.error(error));
    const origin = `http://127.0.0.1:${server.address().port}`;
    let result;
    if (manifest.validation) {
        await page.goto(`${origin}/validation.html`);
        await page.waitForFunction(() => window.__smlibValidationResult, null, {timeout: 120000});
        result = await page.evaluate(() => window.__smlibValidationResult);
    } else {
        await page.goto(`${origin}/`);
        result = await page.evaluate(async () => {
            const api = await import('./index.mjs');
            await Promise.all([api.default(), api.default()]);
            if ('validateSmlib' in api) throw new Error('Validation code shipped in production');
            const backend = new api.SmlibBackend();
            const curve = backend.circle([0,0,0], 2);
            backend.free();
            const value = curve.evaluate(curve.domain()[0]);
            curve.free();
            if (value.length !== 9 || value[0] !== 2) throw new Error('Production loader/geometry failed');
            return {ok:true, production:true};
        });
    }
    console.log(JSON.stringify({browser: await browser.version(), ...result}));
    if (!result.ok) process.exitCode = 1;
} finally {
    if (browser) await browser.close();
    server.close();
}
