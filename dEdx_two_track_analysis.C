
void dEdx_two_track_analysis() {

    const string filenames[] = {
        "./ntuples/TOTEM40.root",
        "./ntuples/TOTEM41.root",
        "./ntuples/TOTEM42.root",
        "./ntuples/TOTEM43.root",
        "./ntuples/TOTEM20.root",
        "./ntuples/TOTEM21.root",
        "./ntuples/TOTEM22.root",
        "./ntuples/TOTEM23.root",
    };

    auto dEdx_hist = new TH2F("dEdx_hist", ";p (GeV);dE/dx",200,0,4,200,0,12);
    int charge_breakers = 0;

    //gStyle->SetOptStat(0);
    gStyle->SetPalette(kCividis);

    for (string filename : filenames) {
        
        cout << "Reading: " << filename << endl;

        TFile *file = TFile::Open(filename.c_str());
        TTree* tree = (TTree*)file->Get("tree");
        TTreeReader Reader(tree);

        TTreeReaderArray<Float_t> trk_p(Reader, "trk_p");
        TTreeReaderArray<Float_t> trk_dedx(Reader, "trk_dedx");
        TTreeReaderValue<Int_t> ntrk(Reader, "ntrk");
        TTreeReaderArray<Int_t> trk_q(Reader, "trk_q");

        while (Reader.Next()) {

            if (trk_p.GetSize() != *ntrk) {
                cout << "ALERT!" << endl;
            }
            if (*ntrk != 4) {
                continue;
            }

            int total_q = 0;
            for (int i=0;i<4;++i){
                Int_t q = trk_q[i];
                total_q += q;
            }

            if (total_q != 0) {
                ++charge_breakers;
                continue;
            }

            for (int i=0;i<4;++i){
                Float_t p = trk_p[i];
                Float_t dedx = trk_dedx[i];
                dEdx_hist->Fill(p, dedx);
            }
        }
    }

    cout << "Charge breakers: " << charge_breakers << endl;

    auto base_canvas = new TCanvas("base_canvas","base_canvas");
    gPad->SetLogz();

    dEdx_hist->Draw("Colz");
}