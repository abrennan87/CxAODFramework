loadrungridjob(const std::string& sampleName)
{
  cout << "Running with ROOT version " << gROOT->GetVersion();
  cout << " (" << gROOT->GetVersionDate() << ")\n";
  cout << "Loading RootCore packages\n";

  gROOT->ProcessLine(".x $ROOTCOREBIN/scripts/load_packages.C");

  cout << "Loading EventLoop grid job\n";

  EL::GridJobLoader gjl;
  gjl.Run(sampleName);
}
