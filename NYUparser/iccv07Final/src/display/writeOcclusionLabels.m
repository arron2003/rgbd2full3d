function writeOcclusionLabels(imdir, segdir, inext, outdir, fn, outext)

for f = 1:numel(fn)
    bn = strtok(fn{f}, '.');
    inname = [segdir '/' bn inext '.mat'];
    load(inname);
    
    outname = [outdir '/' bn outext '.jpg'];
    im = im2double(imread([imdir fn{f}]));
       
    lab = bndinfo.edges.boundaryType;
    
    printOcclusionResult(im, bndinfo, lab, outname, 1);
end