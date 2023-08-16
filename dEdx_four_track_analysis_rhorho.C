
class Particle {
    public:
        //unknown=0, pion=1, rho=2, glueball=3
        int type;
        float p;
        float p_t;
        int charge;
        float eta;
        float phi;
        float dEdx;
        float p_x;
        float p_y;
        float p_z;
        float m;
        float E;

        float dxy;
        float dz;

        Particle(int type_in) {
            type = type_in;   
            m = 0;
        }

        void calculate_3d_momentum() {
            if (p_t==0||eta==0||phi==0) {
                throw invalid_argument("required variables not set");
            }
            p_x = p_t*cos(phi);
            p_y = p_t*sin(phi);
            p_z = p_t*sinh(eta);
        }

        void calculate_energy() {
            if (m == 0) {
                throw invalid_argument("mass of particle not set");
            }
            E = sqrt(pow(p, 2) + pow(m, 2));
        }

        void calculate_total_momentum() {
            p = sqrt(pow(p_x, 2)+pow(p_y, 2)+pow(p_z, 2));
        }

        void calculate_mass() {
            m = sqrt(pow(E, 2) - pow(p, 2));
        }
};

class Event {
    public:
        const float prot_momentum = 6.5e+3;

        Particle** particles;
        int particle_count;

        float ThxR;
        float ThyR;
        float ThxL;
        float ThyL;

        float PR_p[3] = {0,0,prot_momentum};
        float PL_p[3] = {0,0,prot_momentum};


        float diff_p[3] = {0,0,0}; //momentums of the diffractive system: 0=p_x, 1=p_y, 2=p_z

        Event(Particle* particles_in[], int particle_count_in) {
            particles = particles_in;
            particle_count = particle_count_in;
        }

        int calculate_total_charge() {
            float total_charge = 0.;
            for (int i=0; i<particle_count; ++i) {
                Particle* particle = *(particles+i);
                total_charge += particle->charge;
            }
            return total_charge;
        }

        void reconstruct_2_from_4(Particle* origins[2][2], int type) {
            if (particle_count != 4) {
                throw invalid_argument("invalid number of particles");
            }

            Particle* positives[2];
            Particle* negatives[2];

            get_positives_and_negatives(positives, negatives);

            Particle* origin1_c1 = reconstruct_1_from_2(positives[0], negatives[0], type);
            origins[0][0] = origin1_c1;
            Particle* origin2_c1 = reconstruct_1_from_2(positives[1], negatives[1], type);
            origins[0][1] = origin2_c1;
            
            Particle* origin1_c2 = reconstruct_1_from_2(positives[0], negatives[1], type);
            origins[1][0] = origin1_c2;
            Particle* origin2_c2 = reconstruct_1_from_2(positives[1], negatives[0], type);
            origins[1][1] = origin2_c2;
        }

        void get_positives_and_negatives(Particle* positives[2], Particle* negatives[2]) {
            int positive_index = 0;
            int negative_index = 0;
            for (int i=0; i<4; ++i) {
                Particle* current_particle = *(particles+i);
                if (current_particle->charge == 1) {
                    positives[positive_index] = current_particle;
                    positive_index++;
                } else if (current_particle->charge == -1) {
                    negatives[negative_index] = current_particle;
                    negative_index++;
                } else {
                    throw invalid_argument("particle has invalid charge");
                }
            }
        }
    
        Particle* reconstruct_1_from_2(Particle* particle1, Particle* particle2, int type) {
            Particle* origin = new Particle(type);
            origin->E = particle1->E + particle2->E;
            origin->p_x = particle1->p_x + particle2->p_x;
            origin->p_y = particle1->p_y + particle2->p_y;
            origin->p_z = particle1->p_z + particle2->p_z;
            origin->calculate_total_momentum();
            origin->calculate_mass();
            return origin;
        }

        void calculate_proton_momentums() {
            PR_p[0] = ThxR*prot_momentum;
            PL_p[0] = ThxL*prot_momentum;
            PR_p[1] = ThyR*prot_momentum;
            PL_p[1] = ThyL*prot_momentum;
        }

        void calculate_momentum_of_diffractive_system() {
            for (int i=0; i<particle_count; ++i) {
                Particle* current_particle = *(particles+i);
                diff_p[0] += current_particle->p_x;
                diff_p[1] += current_particle->p_y;
                diff_p[2] += current_particle->p_z;
            }
        }

        float get_total_dxy() {
            float total_dxy = 0;
            for (int i=0; i<particle_count; ++i) {
                Particle* current_particle = *(particles+i);
                total_dxy += current_particle->dxy;
            }
            return total_dxy;
        }

        float get_total_dz() {
            float total_dz = 0;
            for (int i=0; i<particle_count; ++i) {
                Particle* current_particle = *(particles+i);
                total_dz += current_particle->dz;
            }
            return total_dz;
        }

        float get_dxy_variance() {
            float average_dxy = get_total_dxy()/particle_count;
            float variance = 0;
            for (int i=0; i<particle_count; ++i) {
                Particle* current_particle = *(particles+i);
                variance += pow(average_dxy-current_particle->dxy, 2);
            }
            return variance;
        }

        float get_dz_variance() {
            float average_dz = get_total_dz()/particle_count;
            float variance = 0;
            for (int i=0; i<particle_count; ++i) {
                Particle* current_particle = *(particles+i);
                variance += pow(average_dz-current_particle->dz, 2);
            }
            return variance;
        }

        float get_squared_total_dxy() {
            float total_squared_dxy = 0;
            for (int i=0; i<particle_count; ++i) {
                Particle* current_particle = *(particles+i);
                total_squared_dxy += pow(current_particle->dxy, 2);
            }
            return total_squared_dxy;
        }

        float get_squared_total_dz() {
            float total_squared_dz = 0;
            for (int i=0; i<particle_count; ++i) {
                Particle* current_particle = *(particles+i);
                total_squared_dz += pow(current_particle->dz, 2);
            }
            return total_squared_dz;
        }

        float get_greatest_dxy() {
            float maximum = 0;
            for (int i=0; i<particle_count; ++i) {
                Particle* current_particle = *(particles+i);
                float dxy = abs(current_particle->dxy);
                if (maximum < dxy) {
                    maximum = dxy;
                }
            }
            return maximum;
        }

        float get_greatest_dz() {
            float maximum = 0;
            for (int i=0; i<particle_count; ++i) {
                Particle* current_particle = *(particles+i);
                float dz = abs(current_particle->dz);
                if (maximum < dz) {
                    maximum = dz;
                }
            }
            return maximum;
        }
};

const float high_k_curve_a = 0.336;
const float high_k_curve_b = 0.14;
const float high_k_curve_c = 3.04;

const float low_k_curve_a = 0.064;
const float low_k_curve_b = 1.33;
const float low_k_curve_c = 3;

const float kaon_vertical_bar = 1.2;
const float kaon_horizontal_bar = 3;

const float relaxed_kaon_vertical_bar = 1.8;

const float pion_vertical_bar = 0.7;
const float pion_horizontal_bar = 1.5;

const float strict_pion_vertical_bar = 0.5;
const float strict_pion_horizontal_bar = 3.5;

const float unknown_vertical_bar = 1.5;

const float two_track_mass_low_limit = 0.755 - 0.100; //740
const float two_track_mass_high_limit = 0.755 + 0.100;

bool is_kaon(float p, float dEdx) {
    if (dEdx < 0.1) {
        return false;
    }
    if (dEdx < kaon_horizontal_bar) {
        return false;
    }
    if (kaon_vertical_bar < p) {
        return false;
    }
    if (dEdx < -1*low_k_curve_a/(p*p)*log(low_k_curve_b*p*p) + low_k_curve_c) {
        return false;
    }
    if (-1*high_k_curve_a/(p*p)*log(high_k_curve_b*p*p) + high_k_curve_c < dEdx) {
        return false;
    }
    return true;
}

bool is_pion(float p, float dEdx) {
    if (dEdx < 0.1) {
        return false;
    }
    if (dEdx < pion_horizontal_bar) {
        return false;
    }
    if (pion_vertical_bar < p) {
        return false;
    }
    if (-1*low_k_curve_a/(p*p)*log(low_k_curve_b*p*p) + low_k_curve_c < dEdx) {
        return false;
    }
    return true;
}

bool strict_is_pion(float p, float dEdx) {
    if (dEdx < 0.1) {
        return false;
    }
    if (dEdx < strict_pion_horizontal_bar) {
        return false;
    }
    if (strict_pion_vertical_bar < p) {
        return false;
    }
    if (-1*low_k_curve_a/(p*p)*log(low_k_curve_b*p*p) + low_k_curve_c < dEdx) {
        return false;
    }
    return true;
}

bool relaxed_is_kaon(float p, float dEdx) {
    if (dEdx < 0.1) {
        return false;
    }
    if (dEdx < kaon_horizontal_bar) {
        return false;
    }
    if (relaxed_kaon_vertical_bar < p) {
        return false;
    }
    if (dEdx < -1*low_k_curve_a/(p*p)*log(low_k_curve_b*p*p) + low_k_curve_c) {
        return false;
    }
    return true;
}

bool is_unknown(float p, float dEdx) {
    if (dEdx < 0.1) {
        return false;
    }
    if (unknown_vertical_bar < p) {
        return true;
    }
    if (kaon_horizontal_bar < dEdx) {
        return false;
    }
    if (p < pion_vertical_bar) {
        return false;
    }
    if (-1*low_k_curve_a/(p*p)*log(low_k_curve_b*p*p) + low_k_curve_c < dEdx) {
        return false;
    }
    return true;
}

bool is_unknown_kaon(float p, float dEdx) {
    bool clear_kaon = true;
    if (dEdx < 0.1) {
        return false;
    }
    if (dEdx < kaon_horizontal_bar) {
        clear_kaon = false;
    }
    if (dEdx < -1*low_k_curve_a/(p*p)*log(low_k_curve_b*p*p) + low_k_curve_c) {
        clear_kaon = false;
    }
    if (-1*high_k_curve_a/(p*p)*log(high_k_curve_b*p*p) + high_k_curve_c < dEdx) {
        clear_kaon = false;
    }
    if (clear_kaon) {
        return true;
    }
    return is_unknown(p, dEdx);
}

bool is_unknown_pion(float p, float dEdx) {
    return (is_pion(p, dEdx) || is_unknown(p, dEdx));
}

void dEdx_four_track_analysis_rhorho() {

    const float pion_mass = 0.13957039;
    const float kaon_mass = 0.493677;

    const string filenames[] = {
        "./ntuples/TOTEM40.root",
        "./ntuples/TOTEM41.root",
        "./ntuples/TOTEM42.root",
        "./ntuples/TOTEM43.root",
        "./ntuples/TOTEM20.root",
        "./ntuples/TOTEM21.root",
        "./ntuples/TOTEM22.root",
        "./ntuples/TOTEM23.root"
    };
    /*
        "./ntuples/TOTEM40.root",
        "./ntuples/TOTEM41.root",
        "./ntuples/TOTEM42.root",
        "./ntuples/TOTEM43.root",
        "./ntuples/TOTEM20.root",
        "./ntuples/TOTEM21.root",
        "./ntuples/TOTEM22.root",
        "./ntuples/TOTEM23.root",
    */

    const float min_dEdx = 5;
    const float max_dEdx = 6;
    const float min_p = 0.33;
    const float max_p = 0.5;

    auto rho_masses = new TH2F("rho_masses", ";m1 (GeV);m2 (GeV)",200,0.2,1.4,200,0.2,1.4);
    auto two_track_mass_1 = new TH2F("two_track_mass_1", ";m1 (GeV);m2 (GeV)",100,0.2,1.4,100,0.2,1.4);
    auto four_track_mass_1 = new TH1F("four_track_mass_1", ";m (GeV)",60,0,3.5);
    auto four_track_mass_2 = new TH1F("four_track_mass_2", ";m (GeV)",60,1.25,2.5);

    auto dEdx_hist = new TH2F("dEdx_hist", ";p (GeV);dE/dx",200,0,4,200,0,12);
    auto kaon_dEdx_hist = new TH2F("kaon_dEdx_hist", ";p (GeV);dE/dx",200,0,4,200,0,12);
    auto pion_dEdx_hist = new TH2F("pion_dEdx_hist", ";p (GeV);dE/dx",200,0,4,200,0,12);
    auto second_dEdx_hist = new TH2F("second_dEdx_hist", ";p (GeV);dE/dx",200,0,4,200,0,12);

    four_track_mass_1->Sumw2();
    four_track_mass_2->Sumw2();

    gStyle->SetOptStat(0);
    gStyle->SetPalette(kCividis);

    for (string filename : filenames) {
        
        cout << "Reading: " << filename << endl;

        TFile *file = TFile::Open(filename.c_str());
        TTree* tree = (TTree*)file->Get("tree");
        TTreeReader Reader(tree);

        TTreeReaderArray<Float_t> trk_dedx(Reader, "trk_dedx");
        TTreeReaderValue<Int_t> ntrk(Reader, "ntrk");

        TTreeReaderArray<Float_t> trk_p(Reader, "trk_p");
        TTreeReaderArray<Float_t> trk_pt(Reader, "trk_pt");
        TTreeReaderArray<Int_t> trk_q(Reader, "trk_q");
        TTreeReaderArray<Float_t> trk_eta(Reader, "trk_eta");
        TTreeReaderArray<Float_t> trk_phi(Reader, "trk_phi");
        while (Reader.Next()) {

            if (trk_p.GetSize() != *ntrk) {
                cout << "ALERT!" << endl;
                throw invalid_argument("ntrk is wrong!");
            }
            if (*ntrk != 4) {
                continue;
            }

            Particle* particles[4];

            for (int i=0;i<4;++i) {
                Particle* particle = new Particle(1);
                particle->p = trk_p[i];
                particle->p_t = trk_pt[i];
                particle->charge = trk_q[i];
                particle->eta = trk_eta[i];
                particle->phi = trk_phi[i];
                particle->dEdx = trk_dedx[i];
                particle->m = pion_mass;
                particle->calculate_3d_momentum();
                particle->calculate_energy();
                particles[i] = particle;
            }
            
            

            Event current_event(particles, 4);
            
            float total_charge = current_event.calculate_total_charge();
            if (total_charge != 0) {
                continue;
            }

            
            Particle* rhos[2][2];
            current_event.reconstruct_2_from_4(rhos, 2);

            float raw_masses[2];
            for (int i=0; i<2; ++i) {
                for (int j=0; j<2; ++j) {
                    Particle* rho = rhos[i][j];
                    raw_masses[j] = rho->m;
                }
                rho_masses->Fill(raw_masses[0], raw_masses[1]);
            }
            
            for (Particle* particle : particles){
                //cout << "hit " << particle->p << ", " << particle->dEdx << endl;
                dEdx_hist->Fill(particle->p, particle->dEdx);
                if (is_kaon(particle->p, particle->dEdx)) {
                    kaon_dEdx_hist->Fill(particle->p, particle->dEdx);
                }
                if (is_pion(particle->p, particle->dEdx)) {
                    pion_dEdx_hist->Fill(particle->p, particle->dEdx);
                }
            }
            
            Particle* positives[2];
            Particle* negatives[2];
            current_event.get_positives_and_negatives(positives, negatives);

            bool valid = false;
            

            if (is_pion(positives[0]->p, positives[0]->dEdx) && is_pion(positives[1]->p, positives[1]->dEdx)) {
                if (is_unknown_pion(negatives[0]->p, negatives[0]->dEdx) && is_unknown_pion(negatives[1]->p, negatives[1]->dEdx)) {
                    valid = true;
                }
            }
            if (is_pion(negatives[0]->p, negatives[0]->dEdx) && is_pion(negatives[1]->p, negatives[1]->dEdx)) {
                if (is_unknown_pion(positives[0]->p, positives[0]->dEdx) && is_unknown_pion(positives[1]->p, positives[1]->dEdx)) {
                    valid = true;
                }
            }
           
            if (valid) {

                for (int i=0;i<4;++i) {
                    particles[i]->m = pion_mass;
                    particles[i]->calculate_energy();
                }

                Particle* rhos[2][2];
                current_event.reconstruct_2_from_4(rhos, 0);

                float masses[2];
                for (int i=0; i<2; ++i) {
                    for (int j=0; j<2; ++j) {
                        Particle* big_kaon = rhos[i][j];
                        masses[j] = big_kaon->m;
                    }
                    two_track_mass_1->Fill(masses[0], masses[1]);

                    if (two_track_mass_low_limit < masses[0] && masses[0] < two_track_mass_high_limit && two_track_mass_low_limit < masses[1] && masses[1] < two_track_mass_high_limit) {
                        Particle* four_track_origin = current_event.reconstruct_1_from_2(rhos[i][0], rhos[i][1], 0);
                        four_track_mass_2->Fill(four_track_origin->m);
                    }

                }

                Particle* four_track_origin = current_event.reconstruct_1_from_2(rhos[0][0], rhos[0][1], 0);
                four_track_mass_1->Fill(four_track_origin->m);

                four_track_origin = current_event.reconstruct_1_from_2(rhos[1][0], rhos[1][1], 0);
                four_track_mass_1->Fill(four_track_origin->m);
            }

            if (strict_is_pion(positives[0]->p, positives[0]->dEdx) && strict_is_pion(positives[1]->p, positives[1]->dEdx)) {
                second_dEdx_hist->Fill(negatives[1]->p, negatives[1]->dEdx);
                second_dEdx_hist->Fill(negatives[0]->p, negatives[0]->dEdx);
            }

            if (strict_is_pion(negatives[0]->p, negatives[0]->dEdx) && strict_is_pion(negatives[1]->p, negatives[1]->dEdx)) {
                second_dEdx_hist->Fill(positives[1]->p, positives[1]->dEdx);
                second_dEdx_hist->Fill(positives[0]->p, positives[0]->dEdx);
            }
        }
    }

    auto rho_canvas = new TCanvas("rho_canvas","rho_canvas");
    rho_masses->Draw("Colz");

    auto two_track_mass_1_canvas = new TCanvas("two_track_mass_1_canvas","two_track_mass_1_canvas");
    two_track_mass_1->Draw("Colz");

    TLine big_K_line_1 = TLine(two_track_mass_low_limit, two_track_mass_low_limit, two_track_mass_high_limit, two_track_mass_low_limit);
    big_K_line_1.SetLineColor(kRed);
    big_K_line_1.SetLineWidth(3);
    big_K_line_1.DrawClone();

    TLine big_K_line_2 = TLine(two_track_mass_high_limit, two_track_mass_low_limit, two_track_mass_high_limit, two_track_mass_high_limit);
    big_K_line_2.SetLineColor(kRed);
    big_K_line_2.SetLineWidth(3);
    big_K_line_2.DrawClone();
    
    TLine big_K_line_3 = TLine(two_track_mass_low_limit, two_track_mass_high_limit, two_track_mass_high_limit, two_track_mass_high_limit);
    big_K_line_3.SetLineColor(kRed);
    big_K_line_3.SetLineWidth(3);
    big_K_line_3.DrawClone();
    
    TLine big_K_line_4 = TLine(two_track_mass_low_limit, two_track_mass_low_limit, two_track_mass_low_limit, two_track_mass_high_limit);
    big_K_line_4.SetLineColor(kRed);
    big_K_line_4.SetLineWidth(3);
    big_K_line_4.DrawClone();


    auto four_track_mass_1_canvas = new TCanvas("four_track_mass_1_canvas","four_track_mass_1_canvas");
    four_track_mass_1->Draw();

    auto four_track_mass_2_canvas = new TCanvas("four_track_mass_2_canvas","four_track_mass_2_canvas");
    four_track_mass_2->Draw();

    

    auto base_canvas = new TCanvas("base_canvas","base_canvas");
    gPad->SetLogz();

    dEdx_hist->Draw("Colz");


    TF1 high_k_curve1("high_k_curve1", "-[0]/(x*x)*log([1]*x*x)+[2]", 0, 1.5);
    high_k_curve1.SetParameters(high_k_curve_a, high_k_curve_b, high_k_curve_c);
    high_k_curve1.DrawCopy("Same");

    TF1 low_k_curve1("low_k_curve1", "-[0]/(x*x)*log([1]*x*x)+[2]", 0, 1.5);
    low_k_curve1.SetParameters(low_k_curve_a, low_k_curve_b, low_k_curve_c);
    low_k_curve1.DrawCopy("Same");


    
    auto K_dEdx_canvas = new TCanvas("K_dEdx_canvas","K_dEdx_canvas");
    gPad->SetLogz();

    kaon_dEdx_hist->Draw("Colz");

    high_k_curve1.DrawCopy("Same");
    low_k_curve1.DrawCopy("Same");


    auto pion_dEdx_canvas = new TCanvas("pion_dEdx_canvas","pion_dEdx_canvas");
    gPad->SetLogz();

    pion_dEdx_hist->Draw("Colz");

    high_k_curve1.DrawCopy("Same");
    low_k_curve1.DrawCopy("Same");


    auto second_dEdx_canvas = new TCanvas("second_dEdx_canvas","second_dEdx_canvas");
    gPad->SetLogz();

    second_dEdx_hist->Draw("Colz");

    high_k_curve1.DrawCopy("Same");
    low_k_curve1.DrawCopy("Same");
}