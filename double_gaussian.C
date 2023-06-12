
void custom_fit_1() {

  auto mycanvas = new TCanvas();

  TF1 gaussian3("gaussian3", "[0]*TMath::Gaus(x,[1],[2])+[3]*TMath::Gaus(x,[4],[5])", -10, 110);
  gaussian3.SetParameters(1,25,10,1,75,10);
  
  TH1F histo("histo","Histogram;X;Y",100,0,100);
  
  histo.Sumw2();
  
  for (int i=0;i<10000;i++) {
    histo.Fill(gaussian3.GetRandom());
  }
  
  
  double total = histo.Integral(1,100);
  
  cout << total << endl;
  
  
  TF1 gaussian_fit_1("gaussian_fit_1", "[0]*TMath::Gaus(x,[1],[2])", 0, 100);
  gaussian_fit_1.SetParameters(1,20,5);
  
  auto fitPtr_1 = histo.Fit(&gaussian_fit_1, "S","",0,50);
  
  
  TF1 gaussian_fit_2("gaussian_fit_2", "[0]*TMath::Gaus(x,[1],[2])", 0, 100);
  gaussian_fit_2.SetParameters(1,20,5);
  
  auto fitPtr_2 = histo.Fit(&gaussian_fit_2, "S","",50,100);
  
  

  //TF1 double_gaussian_fit("double_gaussian_fit", "TMath::Gaus(x,[0],[1])+TMath::Gaus(x,[2],[3])", 0, 100);
  //gaussian_fit.SetParameters(25,10,75,10);
  
  //auto fitPtr_2 = histo.Fit(&double_gaussian_fit, "S");
  
  histo.DrawClone();
  
  //TF1 new_gaussian("gaussian_fit", "[0]*TMath::Gaus(x,[1],[2])", 0, 100);
  
  /*
  for (int i=0;i<3;i++) {
    new_gaussian.SetParameter(i, gaussian_fit.GetParameter(i));
  }
  
  for (int i=0;i<3;i++) {
    cout << new_gaussian.GetParameter(i) << endl;
  } 
  */
 
  gaussian_fit_1.DrawClone("Same");
  gaussian_fit_2.DrawClone("Same");
  
  mycanvas->Print("custom_fit_1.pdf");
  
}
