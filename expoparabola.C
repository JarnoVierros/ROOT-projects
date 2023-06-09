

void expoparabola() {

  auto mycanvas = new TCanvas();
  
  TF1 exponential("exponential", "[0]*TMath::Exp([1]*x)", 0, 100);
  exponential.SetParameters(100, -0.02);
  
  TF1 parabola("parabola", "[0] + [1]*x + [2]*x*x", 40, 60);
  parabola.SetParameters(-2400, 100, -1);
  
  TH1F histo("histo","Histogram;X;Y",100,0,100);
  
  
  for (int i=0;i<10000;i++) {
    histo.Fill(exponential.GetRandom());
  }
  
  for (int i=0;i<1000;i++) {
    histo.Fill(parabola.GetRandom());
  }
  
  TF1 parabola_fit("parabola_fit", "[0] + [1]*x + [2]*x*x", 40, 60);
  parabola_fit.SetParameters(1, 1, -1);
  
  histo.Fit(&parabola_fit, "","",40,60);
  
  string function_string = "[0]*TMath::Exp([1]*x)+((40<x) ? 1 : 0)*((x<60) ? 1 : 0)*(" + to_string(parabola_fit.GetParameter(0)) + "+(" + to_string(parabola_fit.GetParameter(1)) + ")*x+(" + to_string(parabola_fit.GetParameter(2)) + ")*x*x - [0]*TMath::Exp([1]*x))";
  
  TString function_TString = function_string;
  
  cout << "function string: " << function_TString << endl;
  
  TF1 gaussian_fit("gaussian_fit", function_TString, 0, 100);
  gaussian_fit.SetParameters(100, -0.02);
  
  histo.Fit(&gaussian_fit, "","",0,100);
  
  histo.DrawClone();
  gaussian_fit.DrawClone("Same");
  
  
  mycanvas->Print("expoparabola.pdf");
  
  TFile out_file("expoparabola.root"); //,"RECREATE"
  histo.Write();
  
  
  //histo.SaveAs("expoparabola_histogram.root");
  //out_file.Close();
  
  //TCanvas c;
  //histo.Draw();
  //c.SaveAs("histo.pdf");
  
  
}
