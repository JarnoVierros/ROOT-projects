	

void read_data(string filename, double data[]) {
  string line;
  ifstream file(filename);
  int i = 0;
  while (getline (file, line)) { 
    data[i] = stod(line);
    i++;
  }
  file.close();
}

void custom_histogram_1() {

  const int data_length = 20;

  double data[data_length];
  read_data("custom_data_1.txt", data);
  
  for (int i = 0; i < 20; i++) {
    cout << data[i] << endl;
  }
  
  auto histogram = new TH1F("histogram", "Histogram;X;Y", 11, 0, 10);
  
  for (int i = 0; i < data_length; i++) {
    histogram->Fill(data[i]);
  }
  
  
  histogram->Draw();
  
}
