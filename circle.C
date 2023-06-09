

void circle(){

  //TFile ofile("circle.root", "RECREATE");
  TFile ofile("circle.root", "RECREATE");
  
  TTree circle_ntuple("circle_ntuple", "Circle");
  
  float X, Y;
  string quadrant;
  
  circle_ntuple.Branch("X", &X, "X/F");
  circle_ntuple.Branch("Y", &Y, "Y/F");
  circle_ntuple.Branch("quadrant", &quadrant);
  
  TRandom3 rndm;
  
  for (int i=0; i<1000; ++i) {
  
    X = rndm.Uniform(0.,1.);
    Y = sin(acos(X));
  
    int positiveX = rndm.Integer(2);
    if (!positiveX){X = -X;}
  
    int positiveY = rndm.Integer(2);
    if (!positiveY){Y = -Y;}
    
    if (positiveX==1&positiveY==1) {
      quadrant = "topright";
    } else if (positiveX==1&positiveY==0) {
      quadrant = "bottomright";
    } else if (positiveX==0&positiveY==0) {
      quadrant = "bottomleft";
    } else if (positiveX==0&positiveY==1) {
      quadrant = "topleft";
    }
    
    cout << "X: " << X << endl;
    cout << "Y: " << Y << endl;
    cout << "quadrant: " << quadrant << endl;
    
    circle_ntuple.Fill();
  
  }
  
  circle_ntuple.Write();
  ofile.Close();

}
