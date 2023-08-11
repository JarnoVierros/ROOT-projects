
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

void detector_summary_analysis_V2() {

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

    const int detector_count = 2;
    const bool manual_mode = false;
    int skips;

    gStyle->SetPalette(kCividis);

    auto randomgen = TRandom3();
    skips = randomgen.Integer(43703173);
    
    map<int, TH2F*> detectors;
    map<int, float[2]> detector_optimal_eta;

    Int_t max_detector = 0;
    Int_t max_detector_count = 0;

    //auto dEdx = new TH2F("dEdx", ";p;dE/dx", 200,0,2.5,200,0,120);

    const int chosen_ids[] = {344287236, 402664717}; //303054864, 304099348, 470439304, 304189448, 470311846, 402666794
    const int chosen_chips[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
    
    for (string filename : filenames) {
        
        cout << filename << endl;

        TFile *file = TFile::Open(filename.c_str());
        TTree* tree = (TTree*)file->Get("mytree");
        TTreeReader Reader(tree);

        TTreeReaderValue<UInt_t> tree_detector_id(Reader, "tree_detector_id");

        int counter = 0;
        while (Reader.Next()) {
            if (*tree_detector_id != 344287236 && *tree_detector_id != 402664717) {
                continue;
            }
            skips = 0;
            if (counter<skips) { //+50
                //cout << "counter skip "<< counter << "/" << skips << endl;
                ++counter;
                continue;
            }
            cout << "skipped " << skips << endl;
            UInt_t detector_id = *tree_detector_id;
            bool skip = false;
            for (const auto& [key, value] : detectors) {
                if (key == detector_id) {
                    skip = true;
                    break;
                }
            }
            if (skip) {
                cout << "repeat skip" << endl;
                continue;
            }
            TString histogram_name = "histogram_"+to_string(*tree_detector_id);
            detectors[*tree_detector_id] = new TH2F(histogram_name, ";p;dE/dx", 200,0,2.5,200,-30,120);
            cout << "Selected detector: " << *tree_detector_id << endl;

            counter = 0;
            skips = randomgen.Integer(43703173);

            if (detectors.size() >= detector_count) {
                cout << detectors.size() << " detectors selected" << endl;
                break;
            }
            Reader.Restart();
        }
        break;
    }
    
    for (string filename : filenames) {
        
        cout << filename << endl;

        TFile *file = TFile::Open(filename.c_str());
        TTree* tree = (TTree*)file->Get("mytree");
        TTreeReader Reader(tree);

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
            for (const auto& [key, value] : detectors) {
                if (key == *tree_detector_id) {
                    Float_t path_length = 0.1*(*tree_path_length);
                    Float_t p = *tree_momentum;
                    Float_t dE = *tree_deltaE;
                    value->Fill(p, dE/path_length);
                    break;
                }
            }
        }
    }

    cout << "Drawing figures" << endl;

    auto base_canvas = new TCanvas("base_canvas","base_canvas");
    gPad->SetLogz();

    TString filename = "";
    string input;

    for (const auto& [key, value] : detectors) {
        base_canvas->Clear();
        value->Draw("Colz");
        filename = "detector_graphs/"+to_string(key)+".pdf";
        base_canvas->Print(filename);

        if (manual_mode) {
            gPad->Modified();
            gPad->Update();
            base_canvas->Show();

            cout << "Type \"stop\" to close the program. Type anything else to show the projection." << endl;
            input = ""; 
            cin >> input;
            if (input == "stop") {
                cout << "Shutting down." << endl;
                break;
            }

            base_canvas->Clear();

            auto projection = value->ProjectionY("projection");
            projection->Draw();

            gPad->Modified();
            gPad->Update();
            base_canvas->Show();

            cout << "Type \"stop\" to close the program. Type anything else to show the next figure." << endl;
            input = ""; 
            cin >> input;
            if (input == "stop") {
                cout << "Shutting down." << endl;
                break;
            }
        }
    }

    cout << "All done" << endl;

}