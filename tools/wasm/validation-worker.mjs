import load, { validateSmlib } from './index.mjs';
try {
    await load();
    postMessage({ok:true, checks:validateSmlib()});
} catch (error) {
    postMessage({ok:false, error:String(error), stack:error.stack});
}
