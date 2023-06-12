
void ntuple_reader() {

    const string filename = "TOTEM43.root";

    TFile *tree = TFile::Open(filename.c_str());
    TTreeReader Reader("tree", tree);

    TTreeReaderArray<Float_t> trk_eta(Reader, "trk_eta");
    TTreeReaderArray<Float_t> trk_phi(Reader, "trk_phi");

    auto histo = new TH2F("histo", "trk_eta and trk_phi;trk_eta;trk_phi",100,-5,5,100,-5,5);

    while (Reader.Next()) {

        if (trk_eta.GetSize() != 4) {
            continue;
        }

        for (int i=0;i<trk_eta.GetSize();++i) {
            //cout << "trk_eta: " << trk_eta[i] << ", trk_phi: " << trk_phi[i] << endl;
            histo->Fill(trk_eta[i], trk_phi[i]);
        }
        //cout << endl;
    }

    auto main = new TCanvas("Canvas1","Canvas1");
    histo->Draw("Cont4");

    auto projections = new TCanvas("Canvas2","Canvas2");
    projections->Divide(1,2);

    projections->cd(1);
    histo->ProjectionX()->Draw();

    projections->cd(2);
    histo->ProjectionY()->Draw();


}