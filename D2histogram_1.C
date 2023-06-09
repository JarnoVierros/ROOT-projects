

void D2histogram_1() {

  auto histogram = new TCanvas();
  

  TView *view1 = TView::CreateView(1);
  view1->SetRange(0,0,0,100,100,72.765);
  //gStyle->SetPalette(kBird);

  auto histo = new TH2F("histo", "2D Histogram;Exponential;Linear",100,0,100,100,0,100);

  TF1 exp("exp", "100*exp(-0.03*x)", 0, 100);
  TF1 lin("lin", "60-0.5*x", 0, 100);
  
  for (int i = 0; i<100000; i++) {
    histo->Fill(exp.GetRandom(), lin.GetRandom());
  }
  
  
  //mycanvas->SetPhi(-60); // default is 30
  
  histo->DrawClone("Colz");
  
  
  auto projectionY = new TCanvas();

  histo->ProjectionY()->Draw();

  TF1 *lin_fit = new TF1("lin_fit", "[0]+[1]*x",0,100);
  lin_fit->SetParameters(1,-1);
  
  histo->ProjectionY()->Fit(lin_fit);
  

  auto projectionX = new TCanvas();

  histo->ProjectionX()->Draw();

  TF1 *exp_fit = new TF1("exp_fit", "[0]*exp([1]*x)",0,100);
  lin_fit->SetParameters(1,-1);
  
  histo->ProjectionX()->Fit(exp_fit);

}
