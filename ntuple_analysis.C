
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

        Particle(int type_in) {
            type = type_in;   
        }

        void calculate_3d_momentum() {
            if (p_t==0||eta==0||phi==0) {
                throw "required variables not set";
            }
            p_x = p_t*cos(phi);
            p_y = p_t*sin(phi);
            p_z = p_t*sinh(eta);
        };
};

class Event {
    public:
        Particle** particles;
        int particle_count;
        Event(Particle* particles_in[], int particle_count_in) {
            particles = particles_in;
            particle_count = particle_count_in;
        }
};


void ntuple_analysis() {


    //const Particle nul_particle = Particle(0);
    const string filename = "TOTEM43.root";

    TFile *file = TFile::Open(filename.c_str());
    TTree* tree = (TTree*)file->Get("tree");
    TTreeReader Reader(tree);
    //TTreeReader Reader("tree", file);

    TTreeReaderArray<Float_t> trk_p(Reader, "trk_p");
    TTreeReaderArray<Float_t> trk_pt(Reader, "trk_pt");
    TTreeReaderArray<Int_t> trk_q(Reader, "trk_q");
    TTreeReaderArray<Float_t> trk_eta(Reader, "trk_eta");
    TTreeReaderArray<Float_t> trk_phi(Reader, "trk_phi");

    int event_count = tree->GetEntries();
    Event* events[4];

    int event_number = 0;
    while (Reader.Next()) {

        if (trk_p.GetSize() != 4) {
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
            particle->calculate_3d_momentum();
            particles[i] = particle;
        }
        /*
        for (Particle* particle : particles) {
            //cout << "momentum: " << particle.p_t << ", charge: " << particle.charge << endl;
            cout << "p_x: " << particle->p_x << ", p_y: " << particle->p_y << ", p_z: " << particle->p_z << ", Dp: " 
            << particle->p - sqrt(pow(particle->p_x, 2)+pow(particle->p_y, 2)+pow(particle->p_z, 2)) 
            << ", Dp_t: " << particle->p_t - sqrt(pow(particle->p_x, 2)+pow(particle->p_y, 2)) << ", charge: " << particle->charge << endl;
        }
        */
        /*
        cout << *particles << endl;
        cout << particles[0] << endl;
        cout << *(particles + 1) << endl;
        cout << particles[1] << endl;
        cout << particles[2] << endl;
        */
        if (event_number<4) {
            Event current_event(particles, 4);
            events[event_number] = &current_event;
        }


        //cout << endl;
        ++event_number;
    }
    cout << "\n\n\n\nSTART\n";
    
    
    for (Event* event : events) {
        for (int i=0; i<event->particle_count; ++i) {
            Particle* particle = *(event->particles+i);
            cout << particle->p << endl;
        }
        cout << endl;
    }
    
}