word kinematicCloudName("kinematicCloud");
args.optionReadIfPresent("cloudName", kinematicCloudName);

Info<< "Constructing plateletCloud " << kinematicCloudName << endl;
basicPlateletCloud kinematicCloud
(
    kinematicCloudName,
    rho,
    U,
    mu,
    g,
	tau
);
