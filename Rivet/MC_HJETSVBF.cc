// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/NonPromptFinalState.hh"

namespace Rivet {


  /// Higgs + jets study in VBF topology
  class MC_HJETSVBF : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(MC_HJETSVBF);


    /// Book histograms and initialise projections before the run
    void init() {

      // Projections
      NonPromptFinalState fs(Cuts::abseta < 5);
      FastJets fj(fs, FastJets::ANTIKT, 0.4);
      declare(fj, "Jets");
      PromptFinalState pfs(Cuts::abseta < 5);
      PromptFinalState higgses(Cuts::pid == PID::HIGGS && Cuts::abseta < 5);
      declare(pfs, "PromptFS");
      declare(higgses, "Higgses");


      // Histograms, starting with standard m12, dy12 combinations
      /// @todo Add an extra loop, array layer, and suffix element for the jet radius
      for (size_t i = 0; i < SELNAMES.size(); ++i) {
        const string sn = SELNAMES[i];
        book(_c_xs[i], "cross"+(sn.empty() ? "_full" : sn));
        book(_h[i]["deltaphi_jj"], "deltaphi_jj"+sn, 6, 0.0, M_PI);
        book(_h[i]["y_j12_02bin"], "y_j12_02bin"+sn, linspace(20, 0.0, 4.0)+linspace(4, 4.5, 6.5)+linspace(3, 7.0, 10));
        book(_h[i]["pth_largebin"], "pth_largebin"+sn, linspace(4, 0, 100)+linspace(7, 150, 500)+linspace(2, 600, 1000));
        book(_h[i]["pth_finebin"], "pth_finebin"+sn, linspace(19, 0, 95)+linspace(9, 100, 190)+linspace(16, 200, 600));
        book(_h[i]["njets"], "njets"+sn, 5, -0.5, 4.5);
        book(_h[i]["HT"], "HT"+sn, linspace(10, 0, 1000)+linspace(4, 1200, 2000));
        book(_h[i]["pthjj_largebin"], "pthjj_largebin"+sn, linspace(4, 0, 100)+linspace(7, 150, 500)+linspace(2, 600, 1000));
        book(_h[i]["ptjj_largebin"], "ptjj_largebin"+sn, linspace(4, 0, 100)+linspace(7, 150, 500)+linspace(2, 600, 1000));
      }
      // Now the m12 histograms, with various binnings and dy12 cuts
      book(_h_mjj["mjj_STXS"], "mjj_STXS", {0.,100.,200.,350.,700.,1000.,1500.,2000.,2500.,3000.});
      book(_h_mjj["mjj_ATLAS_CONF_029"], "mjj_ATLAS_CONF_029", {0.,160.,500.,1500.});
      book(_h_mjj["mjj_100GeVbin"], "mjj_100GeVbin", 20, 0, 2000);
      book(_h_mjj["mjj_100GeVbin_center"], "mjj_100GeVbin_center", 20, 0, 2000);
      book(_h_mjj["mjj_100GeVbin_middle"], "mjj_100GeVbin_middle", 20, 0, 2000);
      book(_h_mjj["mjj_100GeVbin_forward"], "mjj_100GeVbin_forward", 20, 0, 2000);
      book(_h_mjj["mjj_100GeVbin_ATLAS"], "mjj_100GeVbin_ATLAS", 20, 0, 2000);
      // Now the dy12 histograms, with various m12 cuts
      book(_h_dyjj["deltay_jj"], "deltay_jj", 10, 0, 10);
      book(_h_dyjj["deltay_jj_light"], "deltay_jj_light", 10, 0, 10);
      book(_h_dyjj["deltay_jj_heavy"], "deltay_jj_heavy", 10, 0, 10);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Get Higgs
      const Particles higgses = apply<ParticleFinder>(event, "Higgses").particles();
      if (higgses.empty()) vetoEvent;
      if (higgses.size() > 1) vetoEvent;
      const Particle higgs = higgses[0];
      const double ptH = higgs.pT();

      // Get jets and leading dijet system
      const Jets jets = apply<FastJets>(event, "Jets").jetsByPt(Cuts::pT > 30*GeV && Cuts::absrap < 4.4);
      if (jets.size() < 2) vetoEvent;
      const int njets = jets.size();
      const double htj = sum(jets, Kin::pT, 0.0);
      const double ht = htj + higgs.pT();

      // Get leading dijet system
      const FourMomentum dijet = jets[0].mom() + jets[1].mom();
      const double m12 = dijet.mass();
      const double pt12 = dijet.pT();
      const double y12 = dijet.rap();
      const double dy12 = abs(jets[0].rap() - jets[1].rap());
      const double dphi12 = deltaPhi(jets[0], jets[1]);

      // Get leading dijet + H
      const FourMomentum jjH = dijet + higgs.mom();
      const double ptjjH = jjH.pT();
      const double dphijjH = deltaPhi(higgs, dijet); //min(deltaPhi(higgs, jets[0]), deltaPhi(higgs, jets[1]));

      // Fill standard cut-combination histograms
      for (size_t i = 0; i < SELNAMES.size(); ++i) {
        if (m12 < M12CUTS[i][0]*GeV) continue;
        if (m12 >= M12CUTS[i][1]*GeV) continue;
        if (dy12 < DY12CUTS[i][0]) continue;
        if (dy12 >= DY12CUTS[i][1]) continue;
        if (dphijjH < DPHIJHCUTS[i]) continue;
        _c_xs[i]->fill();
        _h[i]["deltaphi_jj"]->fill(dphi12);
        _h[i]["y_j12_02bin"]->fill(y12);
        _h[i]["pth_largebin"]->fill(ptH/GeV);
        _h[i]["pth_finebin"]->fill(ptH/GeV);
        _h[i]["njets"]->fill(njets);
        _h[i]["HT"]->fill(ht/GeV);
        _h[i]["pthjj_largebin"]->fill(ptjjH/GeV);
        _h[i]["ptjj_largebin"]->fill(pt12/GeV);
      }

      // Fill dijet mass histograms
      _h_mjj["mjj_STXS"]->fill(m12/GeV);
      _h_mjj["mjj_ATLAS_CONF_029"]->fill(m12/GeV);
      _h_mjj["mjj_100GeVbin"]->fill(m12/GeV);
      if (dy12 < 2)
        _h_mjj["mjj_100GeVbin_center"]->fill(m12/GeV);
      if (dy12 >= 2 && dy12 < 4)
        _h_mjj["mjj_100GeVbin_middle"]->fill(m12/GeV);
      if (dy12 > 4)
        _h_mjj["mjj_100GeVbin_forward"]->fill(m12/GeV);
      if (dy12 > 3 && m12 > 400*GeV && dphijjH > 2.8)
        _h_mjj["mjj_100GeVbin_ATLAS"]->fill(m12/GeV);

      // Now the dy12 histograms, with various m12 cuts
      _h_dyjj["deltay_jj"]->fill(dy12);
      if (m12 < 350*GeV)
        _h_dyjj["deltay_jj_light"]->fill(dy12);
      else // if (m12 < 350*GeV)
        _h_dyjj["deltay_jj_heavy"]->fill(dy12);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for (size_t i = 0; i < SELNAMES.size(); ++i) {
        scale(_c_xs[i], crossSection()/femtobarn/sumOfWeights());
        for (auto kv : _h[i])
          scale(kv.second, crossSection()/femtobarn/sumOfWeights());
      }
      for (auto kv : _h_mjj)
        scale(kv.second, crossSection()/femtobarn/sumOfWeights());
      for (auto kv : _h_dyjj)
        scale(kv.second, crossSection()/femtobarn/sumOfWeights());
    }


    /// @name Histograms
    //@{
    map<string, Histo1DPtr> _h[8];
    CounterPtr _c_xs[8];
    map<string,Histo1DPtr> _h_mjj, _h_dyjj;
    //@}


    /// Cut values for standard histogram sets (other than m12 and dy12)
    //@{
    static const vector<string> SELNAMES;
    static const vector<doubles> M12CUTS;
    static const vector<doubles> DY12CUTS;
    static const vector<double> DPHIJHCUTS;
    //@}

  };


  // Static const initializers
  const vector<string>  MC_HJETSVBF::SELNAMES = {"", "_light_center", "_heavy_center", "_light_middle", "_heavy_middle", "_light_forward", "_heavy_forward", "_ATLAS"};
  const vector<doubles> MC_HJETSVBF::M12CUTS = {{0., HUGE_VAL}, {0., 350.}, {350., HUGE_VAL}, {0., 350.}, {350., HUGE_VAL}, {0., 350.}, {350., HUGE_VAL}, {400., HUGE_VAL}};
  const vector<doubles> MC_HJETSVBF::DY12CUTS = {{0., HUGE_VAL}, {0., 2.}, {2., 4.}, {4., HUGE_VAL}, {0., 2.}, {2., 4.}, {4., HUGE_VAL}, {3., HUGE_VAL}};
  const vector<double>  MC_HJETSVBF::DPHIJHCUTS = {0., 0., 0., 0., 0., 0., 0., 2.8};


  DECLARE_RIVET_PLUGIN(MC_HJETSVBF);

}
