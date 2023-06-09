

void circle_reader(){
  
  
  auto circle_histogram = new TH2F("circle_histogram", "Circle;X;Y",100,-2,2,100,0-2,2);
  
  TFile *myFile = TFile::Open("circle.root");
  TTreeReader myReader("circle_ntuple", myFile);
  
  TTreeReaderValue<Float_t> X(myReader, "X");
  TTreeReaderValue<Float_t> Y(myReader, "Y");
  TTreeReaderValue<string> quadrant(myReader, "quadrant");
  
  while (myReader.Next()) {
    circle_histogram->Fill(*X, *Y);
    cout << "X: " <<  *X << endl;
    cout << "Y: " << *Y << endl;
    cout << "quadrant: " << *quadrant << endl;
    cout << endl;
  }
  
  circle_histogram->Draw();
  
}
