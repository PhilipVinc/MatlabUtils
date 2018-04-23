function output = ReadAnalyzeChunk( obj , i, storeFlag)

	if (obj.chunkData{i} ~= [])
        nLoaded = obj.ReadAllChunk(i);
    end

	fprintf( ['\tAnalizing...']);
   	res = obj.AverageExtractData(obj.chunkData{i}, obj.params);
   	
   	fprintf( ['\t Done!\n\tSaving Data...']);
   	obj.params = res.params;
   	params = obj.params;
   	params.n_traj = obj.chunkTrajN(i);
   	averaged = res.ave;
   	quantities = res.quan;
   	if i==1
   	    save(aveCnkPath, 'averaged', 'params');
   	else
   	    clear params;
   	    params.n_traj = obj.chunkTrajN(i);
   	    save(aveCnkPath, 'averaged', 'params');
   	end            
   	% This fix is needed because thtop has a future date.
   	obj.FixEditDate(aveCnkPath, cIPath);
   	quanCnkPath = obj.QuantitiesChunkPath(chunkId);
   	save(quanCnkPath, 'params', 'quantities');
   	obj.FixEditDate(quanCnkPath, cIPath);
   	fprintf( ['\t Done!\n']);
   	clear params; clear averaged; clear quantities;        

   	if ~storeFlag
   		obj.chunkData{i} = [];
   	end

end