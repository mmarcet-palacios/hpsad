const { chromium } = require('playwright');
(async()=>{const b=await chromium.launch();const p=await b.newPage({viewport:{width:1040,height:1200},deviceScaleFactor:3});
for (const [src,out] of process.argv.slice(2).reduce((a,x,i,arr)=>i%2?a:a.concat([[x,arr[i+1]]]),[])){
 await p.goto('file://'+src);await p.waitForTimeout(1500);
 const el=await p.$('svg'); await el.screenshot({path:out});}
await b.close();})();
