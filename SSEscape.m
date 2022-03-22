function [BGDepth,BGSigma,EscapeFeatures] = SSEscape(PG,BGSeg,BGSigma,BGPhi,BGDepth,EscapeFeatures)
BGFeature = [PG.get('edgeNum',BGSeg(1,1));PG.get('edgeNum',BGSeg(1,1))
featurepair = BGFeature;
featurepair(1) = 0;
[ppv,sigma] = depthperfeature2(PG,BGSigma,BGPhi,featurepair);

UniqueEscapeFeatures = unique(EscapeFeatures);
end