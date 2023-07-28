
#include "CMS_lumi.C"

string create_interval(float start, float stop) {
    return "((("+to_string(start)+"<x) ? 1 : 0) - (("+to_string(stop)+"<x) ? 1 : 0))";
}

int enforce_interval(float x, float start, float stop) {
    return ((start<x ? 1 : 0) - (stop<x ? 1 : 0));
}

string get_par(TF1 fit, int par) {
    return to_string(fit.GetParameter(par));
}

class Particle {
    public:
        //pion=1, rho=2
        int type;
        float p;
        float p_t;
        int charge;
        float eta;
        float phi;
        float p_x;
        float p_y;
        float p_z;
        float m;
        float E;

        float dxy;
        float dz;

        Particle(int type_in) {
            type = type_in;   
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
        const float prot_momentum = 6.5e+6;

        Particle** particles;
        int particle_count;

        float ThxR;
        float ThyR;
        float ThxL;
        float ThyL;

        float PR_p[3] = {0,0,prot_momentum};
        float PL_p[3] = {0,0,prot_momentum};


        float ref_p[3] = {0,0,0}; //momentums of the refractive system: 0=p_x, 1=p_y, 2=p_z

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

        void reconstruct_2_rhos(Particle* rhos[2][2]) {
            if (particle_count != 4) {
                throw invalid_argument("invalid number of particles");
            }

            Particle* positives[2];
            int positive_index = 0;
            Particle* negatives[2];
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

            Particle* rho1_c1 = reconstruct_rho(positives[0], negatives[0]);
            rhos[0][0] = rho1_c1;
            Particle* rho2_c1 = reconstruct_rho(positives[1], negatives[1]);
            rhos[0][1] = rho2_c1;
            
            Particle* rho1_c2 = reconstruct_rho(positives[0], negatives[1]);
            rhos[1][0] = rho1_c2;
            Particle* rho2_c2 = reconstruct_rho(positives[1], negatives[0]);
            rhos[1][1] = rho2_c2;
        }
    
        Particle* reconstruct_rho(Particle* positive, Particle* negative) {
            Particle* rho = new Particle(2);
            rho->E = positive->E + negative->E;
            rho->p_x = positive->p_x + negative->p_x;
            rho->p_y = positive->p_y + negative->p_y;
            rho->p_z = positive->p_z + negative->p_z;
            rho->calculate_total_momentum();
            rho->calculate_mass();
            return rho;
        }

        void calculate_proton_momentums() {
            PR_p[0] = ThxR*prot_momentum;
            PL_p[0] = ThxL*prot_momentum;
            PR_p[1] = ThyR*prot_momentum;
            PL_p[1] = ThyL*prot_momentum;
        }

        void calculate_momentum_of_refractive_system() {
            for (int i=0; i<particle_count; ++i) {
                Particle* current_particle = *(particles+i);
                ref_p[0] += current_particle->p_x;
                ref_p[1] += current_particle->p_y;
                ref_p[2] += current_particle->p_z;
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
};

void refractive_momentum_analysis() {

    //const string filenames[4] = {"./ntuples/TOTEM20.root", "./ntuples/TOTEM21.root", "./ntuples/TOTEM22.root", "./ntuples/TOTEM23.root"};
    //const string filenames[4] = {"./ntuples/TOTEM40.root", "./ntuples/TOTEM41.root", "./ntuples/TOTEM42.root", "./ntuples/TOTEM43_old.root"};
    const string filenames[1] = {"./ntuples/TOTEM43.root"};
    
    const float pion_mass = 139.57039;

    auto rho_masses = new TH2F("rho_masses", "Masses of potential rho particles;m1/MeV;m2/MeV",200,200,1400,200,200,1400);

    auto refractive_momentums = new TH2F("refractive_momentums", ";px (MeV);py (MeV)",200,-1400,1400,200,-1400,1400);

    gStyle->SetOptStat(0);
    gStyle->SetPalette(kCividis);
    rho_masses->Sumw2();

    for (string filename : filenames) {

        TFile *file = TFile::Open(filename.c_str());
        TTree* tree = (TTree*)file->Get("tree");
        TTreeReader Reader(tree);
        //TTreeReader Reader("tree", file);

        TTreeReaderArray<Float_t> trk_p(Reader, "trk_p");
        TTreeReaderArray<Float_t> trk_pt(Reader, "trk_pt");
        TTreeReaderArray<Int_t> trk_q(Reader, "trk_q");
        TTreeReaderArray<Float_t> trk_eta(Reader, "trk_eta");
        TTreeReaderArray<Float_t> trk_phi(Reader, "trk_phi");
        

        while (Reader.Next()) {

            if (trk_p.GetSize() != 4) {
                continue;
            }

            Particle* particles[4];

            for (int i=0;i<4;++i) {
                Particle* particle = new Particle(1);
                particle->p = 1000*trk_p[i];
                particle->p_t = 1000*trk_pt[i];
                particle->charge = trk_q[i];
                particle->eta = trk_eta[i];
                particle->phi = trk_phi[i];
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
            current_event.reconstruct_2_rhos(rhos);

            current_event.calculate_momentum_of_refractive_system();
            
            refractive_momentums->Fill(current_event.ref_p[0], current_event.ref_p[1]);
            

            float masses[2];
            for (int i=0; i<2; ++i) {
                for (int j=0; j<2; ++j) {
                    Particle* rho = rhos[i][j];
                    masses[j] = rho->m;
                    //cout << "m: " << rho->m << endl;
                }
                rho_masses->Fill(masses[0], masses[1]);
            }
        }
    }


    auto m1_vs_m2 = new TCanvas("Canvas1","Canvas1");
    rho_masses->Draw("Colz");


    auto projection_canvas = new TCanvas("Canvas2","Canvas2");
    //auto projection = rho_masses->ProjectionX("X_projection", (rho_mass-allowed_rho_2_mass_difference-200)/6, (rho_mass+allowed_rho_2_mass_difference-200)/6);
    auto projection = rho_masses->ProjectionX("X_projection");
    projection->Draw();


    auto refractive_momentums_canvas = new TCanvas("Canvas3","Canvas3");
    refractive_momentums->Draw("Colz");

    CMS_lumi( refractive_momentums_canvas, 17, 11 );

}