
class Detection {
    public:
        UInt_t detector_id;
        Int_t chip_id;
        Float_t momentum;
        Float_t pt;
        Float_t path_length;
        Float_t deltaE;
        Float_t eta;
        Float_t kaon;
        Float_t clean_RP_conf;
        Detection() {  
        }
};

class Track {
    public:
        vector<Detection*> detections;
        Float_t p_t;
        Float_t total_deltaE = 0;
        Float_t total_path_length = 0;
        Float_t total_p = 0;
        Float_t average_p;
        Float_t dEdx;

        Track() {  
        }

        void reset() {
            detections.clear();
        }

        void calculate_all() {
            calculate_total_deltaE();
            calculate_total_path_length();
            calculate_total_p();
            calculate_average_p();
            calculate_dEdx();
        }

        void calculate_total_deltaE() {
            Float_t total = 0;
            for (Detection* detection : detections) {
                total += detection->deltaE;
            }
            total_deltaE = total;
        }

        void calculate_total_path_length() {
            Float_t total = 0;
            for (Detection* detection : detections) {
                total += detection->path_length;
            }
            total_path_length = total;
        }

        void calculate_total_p() {
            Float_t total = 0;
            for (Detection* detection : detections) {
                total += detection->momentum;
            }
            total_p = total;
        }

        void calculate_average_p() {
            if (total_p == 0) {
                calculate_total_p();
            }
            average_p = total_p/detections.size();
        }

        void calculate_dEdx(){
            if (total_deltaE == 0) {
                calculate_total_deltaE();
            }
            if (total_path_length == 0) {
                calculate_total_path_length();
            }
            dEdx = total_deltaE/total_path_length;
        }
};

string set_length(float value, int length) {
    string value_string = to_string(value);
    int value_string_length = value_string.length();
    int difference = length - value_string_length;
    string final_string = "";
    if (difference > 0) {
        for (int i=0; i<difference; ++i) {
            final_string += " ";
        }
        final_string += value_string;
    } else if (difference < 0) {
        for (int i=0; i<value_string_length+difference; i++) {
            final_string += value_string[i];
        }
    } else {
        final_string += value_string;
    }
    return final_string;
}

string set_length(int value, int length) {
    string value_string = to_string(value);
    int difference = length - value_string.length();
    if (difference < 0) {
        throw invalid_argument("length is too short or value is too long"); 
    }
    string final_string = "";
    for (int i=0; i<difference; ++i) {
        final_string += "0";
    }
    final_string += value_string;
    return final_string;
}

void track_figure_summary_analysis() {

    const string filenames[] = {"./ntuples/summary_00.root"};//"./ntuples/summary_01.root"

    auto dEdx = new TH2F("dEdx", ";p;dEdx",200,0,2.5,200,0,120);
    auto path_vs_dE = new TH2F("path_vs_dE", ";path length;delta E",200,0,0.01,200,0,0.5);
    auto track_dEdx = new TH1F("track_dEdx", ";dE/dx;detections",200,0,100);
    
    //map<int, int> detectors;

    Float_t maximum = 0;
    int counter = 0;

    for (string filename : filenames) {
        
        cout << filename << endl;

        TFile *file = TFile::Open(filename.c_str());
        TTree* tree = (TTree*)file->Get("mytree");
        TTreeReader Reader(tree);

        //gStyle->SetOptStat(0);
        //gStyle->SetPalette(kCividis);


        TTreeReaderValue<UInt_t> tree_detector_id(Reader, "tree_detector_id");
        TTreeReaderValue<Int_t> tree_chip_id(Reader, "tree_chip_id");
        TTreeReaderValue<Float_t> tree_momentum(Reader, "tree_momentum");
        TTreeReaderValue<Float_t> tree_pt(Reader, "tree_pt");
        TTreeReaderValue<Float_t> tree_path_length(Reader, "tree_path_length");
        TTreeReaderValue<Float_t> tree_deltaE(Reader, "tree_deltaE");
        TTreeReaderValue<Float_t> tree_eta(Reader, "tree_eta");
        TTreeReaderValue<Bool_t> tree_kaon(Reader, "tree_kaon");
        TTreeReaderValue<Bool_t> tree_clean_RP_conf(Reader, "tree_clean_RP_conf");

        float previous_tree_pt = 0;

        Track* track = new Track();


        while (Reader.Next()) {
            if (*tree_pt != previous_tree_pt && previous_tree_pt != 0) {

                track->calculate_all();

                if (-1 < counter) {
                    for (Detection* detection : track->detections) {
                        path_vs_dE->Fill(detection->path_length, detection->deltaE);
                        track_dEdx->Fill((detection->deltaE)/(detection->path_length));
                    }

                    auto base_canvas = new TCanvas("base_canvas","base_canvas");
                    path_vs_dE->Draw("Colz");

                    auto projection = path_vs_dE->ProjectionY("projection");
                    auto projection_canvas = new TCanvas("projection_canvas","projection_canvas");
                    projection->Draw();

                    auto track_dEdx_canvas = new TCanvas("track_dEdx_canvas","track_dEdx_canvas");
                    track_dEdx->Draw("Colz");

                    unique_ptr<TFile> rootFile(TFile::Open("track_dEdx.root", "RECREATE"));
                    rootFile->WriteObject(&track_dEdx, "track_dEdx");

                    break;
                }
                ++counter;

                dEdx->Fill(track->average_p, track->dEdx);
                track->reset();
                
            }
            
            previous_tree_pt = *tree_pt;

            
            

            Detection* detection = new Detection();
            detection->detector_id = *tree_detector_id;
            detection->chip_id = *tree_chip_id;
            detection->momentum = *tree_momentum;
            detection->pt = *tree_pt;
            detection->path_length = 0.1*(*tree_path_length);
            detection->deltaE = *tree_deltaE;
            detection->eta = *tree_eta;
            detection->kaon = *tree_kaon;
            detection->clean_RP_conf = *tree_clean_RP_conf;


            track->detections.push_back(detection);
            

        }
    }



    //auto base_canvas = new TCanvas("base_canvas","base_canvas");

    //dEdx->Draw("Colz");


}