function output = FixEditDate( obj, aveCnkPath, dataPath )

	aveObj = dir(aveCnkPath);
    try
        dateAve = datetime(aveObj(1).date, 'Locale',  get(0, 'Language'));
    catch
        dateAve = datetime(aveObj(1).date, 'Locale',  'en');
    end
    dataObj = dir(dataPath);
    try
        dateData = datetime(dataObj(1).date, 'Locale',  get(0, 'Language'));
    catch
        dateData = datetime(dataObj(1).date, 'Locale',  'en');
    end
    % Means we  created the data in the future? nopey dopey.
    if (dateData > dateAve)
    	system(['touch -r ', dataPath, ' ', aveCnkPath]);
    end
end