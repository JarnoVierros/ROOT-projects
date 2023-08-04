
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

void detector_summary_analysis() {

    const string filenames[] = {
        "./ntuples/summary_00.root",
        "./ntuples/summary_01.root",
        "./ntuples/summary_02.root",
        "./ntuples/summary_03.root",
        "./ntuples/summary_04.root",
        "./ntuples/summary_05.root",
        "./ntuples/summary_06.root",
        "./ntuples/summary_07.root",
        "./ntuples/summary_08.root",
        "./ntuples/summary_09.root",
        "./ntuples/summary_10.root",
        "./ntuples/summary_11.root",
        "./ntuples/summary_12.root",
        "./ntuples/summary_13.root",
        "./ntuples/summary_14.root",
    };
    
    vector<int> detectors;
    map<int, float[2]> detector_optimal_eta;

    Int_t max_detector = 0;
    Int_t max_detector_count = 0;

    auto dEdx = new TH2F("dEdx", ";p;dE/dx", 200,0,2.5,200,0,120);

    //const int chosen_ids[] = {470439304}; //303054864, 304099348, 470439304, 304189448, 470311846, 402666794
    const int chosen_chips[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};

    int detectors_number = 0;

    for (string filename : filenames) {
        
        cout << filename << endl;

        TFile *file = TFile::Open(filename.c_str());
        TTree* tree = (TTree*)file->Get("mytree");
        TTreeReader Reader(tree);

        //gStyle->SetOptStat(0);
        //gStyle->SetPalette(kCividis);

        TTreeReaderValue<UInt_t> tree_detector_id(Reader, "tree_detector_id");

        int counter = 0;
        while (Reader.Next()) {
            /*
            if (counter<10) {
                ++counter;
                continue;
            }
            */
            UInt_t detector_id = *tree_detector_id;
            bool skip = false;
            for (int detector : detectors) {
                if (detector == detector_id) {
                    skip = true;
                    break;
                }
            }
            if (skip) {
                continue;
            }
            
            ++detectors_number;
            
            continue;

            detectors.push_back(*tree_detector_id);

            if (detectors.size() > 10) {
                break;
            }
        }
        break;
    }
    cout << "Detectors: " << detectors_number << endl;
    return;

    auto base_canvas = new TCanvas("base_canvas","base_canvas");

    for (UInt_t detector : detectors) {

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


            while (Reader.Next()) {

                if (*tree_detector_id != detector) {
                    continue;
                }   
         
                /*
                bool chosen = false;
                for (int chosen_id : chosen_ids){
                    if (detector_id == chosen_id) {
                        chosen = true;
                    }
                }
                if (!chosen) {
                    continue;
                }

                chosen = false;
                for (int chosen_chip : chosen_chips){
                    if (chip_id == chosen_chip) {
                        chosen = true;
                    }
                }
                if (!chosen) {
                    continue;
                }
                */

                Float_t path_length = 0.1*(*tree_path_length);
                Float_t p = *tree_momentum;
                Float_t dE = *tree_deltaE;

                dEdx->Fill(p, dE/path_length);
            }
        }

        /*
        cout << "Type \"stop\" to close the program. Type anything else to show the next figure." << endl;
        string input = ""; 
        cin >> input;

        if (input == "stop") {
            cout << "Shutting down." << endl;
            break;
        }
        */

        cout << "Detector: " << detector << endl;

        base_canvas->Clear();

        dEdx->Draw("Colz");
        gPad->Modified();
        gPad->Update();
        base_canvas->Show();

        TString filename = "detector_graphs/"+to_string(detector)+".pdf";
        base_canvas->Print(filename);

        dEdx->Reset();

        gSystem->ProcessEvents();
    
        cout << "Preparing next figure." << endl;

    }
//303075348
    //path_length_dev_dE->Draw("Colz");
    //detector_eta_hist->Draw();
    //etadE->Draw("Colz");
    //dEdx->Draw("Colz");
    //detector_ids->Draw();

}