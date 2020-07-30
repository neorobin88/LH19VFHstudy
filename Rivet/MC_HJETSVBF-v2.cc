// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/NonPromptFinalState.hh"
using namespace std;
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
      for (size_t i = 0; i < DRS.size(); ++i) {
	FastJets fj(fs, FastJets::ANTIKT, DRS[i]);
	std::string dr = "_xxx";
	sprintf(&dr[1], "%0.1f", DRS[i]);
	declare(fj, "Jets"+dr);
      }
      PromptFinalState pfs(Cuts::abseta < 5);
      PromptFinalState higgses(Cuts::pid == PID::HIGGS && Cuts::abseta < 5);
      declare(pfs, "PromptFS");
      declare(higgses, "Higgses");


      // Histograms, starting with standard m12, dy12 combinations
      /// @todo Add an extra loop, array layer, and suffix element for the jet radius
      for (size_t j = 0; j < DRS.size(); ++j) {
	std::string dr = "_xxx";
	sprintf(&dr[1], "%0.1f", DRS[j]);
	for (size_t i = 0; i < SELNAMES.size(); ++i) {
	  const string sn = dr+SELNAMES[i];
	  book(_c_xs[j][i], "cross"+sn);
	  book(_h[j][i]["deltaphi_jj"], "deltaphi_jj"+sn, 6, 0.0, M_PI);
	  book(_h[j][i]["y_j12_02bin"], "y_j12_02bin"+sn, linspace(20, 0.0, 4.0)+linspace(4, 4.5, 6.5)+linspace(3, 7.0, 10));
	  book(_h[j][i]["pth_largebin"], "pth_largebin"+sn, linspace(4, 0, 100)+linspace(7, 150, 500)+linspace(2, 600, 1000));
	  book(_h[j][i]["pth_finebin"], "pth_finebin"+sn, linspace(19, 0, 95)+linspace(9, 100, 190)+linspace(16, 200, 600));
	  book(_h[j][i]["njets"], "njets"+sn, 5, -0.5, 4.5);
	  book(_h[j][i]["HT"], "HT"+sn, linspace(10, 0, 1000)+linspace(4, 1200, 2000));
	  book(_h[j][i]["pthjj_largebin"], "pthjj_largebin"+sn, linspace(4, 0, 100)+linspace(7, 150, 500)+linspace(2, 600, 1000));
	  book(_h[j][i]["logpthjj_largebin"], "logpthjj_largebin"+sn, linspace(100, 0., 4.));
	  book(_h[j][i]["ptjj_largebin"], "ptjj_largebin"+sn, linspace(4, 0, 100)+linspace(7, 150, 500)+linspace(2, 600, 1000));
	  // for (size_t c = 0; c < HPTCUTS.size(); ++c) {
	  //   std::string pt = "_xxx";
	  //   sprintf(&pt[1], "%0.0f", HPTCUTS[c]);
	  //   book(_h[j][i]["deltay_jj"+pt], "deltay_jj"+sn+pt, 100, 0, 10);
	  //   book(_h[j][i]["deltaR_jj"+pt], "deltaR_jj"+sn+pt, 100, 0, 10);
	  //   book(_h[j][i]["m_jj"+pt], "m_jj"+sn+pt, 40, 0, 2000);
	  // }
	}
	// Now the m12 histograms, with various binnings and dy12 cuts
	book(_h_mjj[j]["mjj_STXS"], "mjj_STXS"+dr, {0.,100.,200.,350.,700.,1000.,1500.,2000.,2500.,3000.});
	book(_h_mjj[j]["mjj_ATLAS_CONF_029"], "mjj_ATLAS_CONF_029"+dr, {0.,160.,500.,1500.});
	book(_h_mjj[j]["mjj_100GeVbin"], "mjj_100GeVbin"+dr, 20, 0, 2000);
	book(_h_mjj[j]["mjj_100GeVbin_center"], "mjj_100GeVbin_center"+dr, 20, 0, 2000);
	book(_h_mjj[j]["mjj_100GeVbin_middle"], "mjj_100GeVbin_middle"+dr, 20, 0, 2000);
	book(_h_mjj[j]["mjj_100GeVbin_forward"], "mjj_100GeVbin_forward"+dr, 20, 0, 2000);
	book(_h_mjj[j]["mjj_100GeVbin_ATLAS"], "mjj_100GeVbin_ATLAS"+dr, 20, 0, 2000);
	// // Now the dy12 histograms, with various m12 cuts
	// book(_h_dyjj[j]["deltay_jj"], "deltay_jj"+dr, 10, 0, 10);
	// book(_h_dyjj[j]["deltay_jj_light"], "deltay_jj_light"+dr, 10, 0, 10);
	// book(_h_dyjj[j]["deltay_jj_heavy"], "deltay_jj_heavy"+dr, 10, 0, 10);

	for (size_t c = 0; c < HPTCUTS.size(); ++c) {
	  std::string pt = std::to_string(int(HPTCUTS[c]));
	  pt = "_pt"+pt;
	  book(_hh[j]["deltay_jj"+pt], "deltay_jj"+dr+pt, 100, 0, 10);
	  book(_hh[j]["deltaR_jj"+pt], "deltaR_jj"+dr+pt, 100, 0, 10);
	  book(_hh[j]["deltaphi_jj"+pt], "deltaphi_jj"+dr+pt, 24, 0, M_PI);
	  book(_hh[j]["y_h"+pt], "y_h"+dr+pt, 70, -3.5, 3.5);
	  book(_hh[j]["y_j1"+pt], "y_h"+dr+pt, 70, -3.5, 3.5);
	  book(_hh[j]["y_j2"+pt], "y_h"+dr+pt, 70, -3.5, 3.5);
	  book(_hh[j]["y_j3"+pt], "y_h"+dr+pt, 70, -3.5, 3.5);
	  book(_hh[j]["m_jj"+pt], "m_jj"+dr+pt, 80, 0, 4000);
	}

	for (size_t c = 0; c < HPTCUTSV2.size(); ++c) {
	  std::string pt = std::to_string(int(HPTCUTSV2[c]));
	  pt = "_pt"+pt;
	  for (size_t c1 = 0; c1 < M12MAXCUTS.size(); ++c1) {
	    std::string mpt = std::to_string(int(M12MAXCUTS[c1]));
	    mpt = "_m"+mpt+pt;
	    book(_hh[j]["deltay_jj"+mpt], "deltay_jj"+dr+mpt, 100, 0, 10);
	    book(_hh[j]["deltaR_jj"+mpt], "deltaR_jj"+dr+mpt, 100, 0, 10);
	    book(_hh[j]["deltaphi_jj"+mpt], "deltaphi_jj"+dr+mpt, 24, 0, M_PI);
	    book(_hh[j]["y_h"+mpt], "y_h"+dr+mpt, 70, -3.5, 3.5);
	    book(_hh[j]["y_j1"+mpt], "y_h"+dr+mpt, 70, -3.5, 3.5);
	    book(_hh[j]["y_j2"+mpt], "y_h"+dr+mpt, 70, -3.5, 3.5);
	    book(_hh[j]["y_j3"+mpt], "y_h"+dr+mpt, 70, -3.5, 3.5);
	  }
	}

	for (size_t c = 0; c < HPTCUTSV2.size(); ++c) {
	  std::string pt = std::to_string(int(HPTCUTSV2[c]));
	  pt = "_pt"+pt;
	  for (size_t c1 = 0; c1 < DR12CUTS.size(); ++c1) {
	    std::string drpt =  "_drxxx";
	    sprintf(&drpt[3], "%0.1f", DR12CUTS[c1]);
	    drpt = drpt+pt;
	    book(_hh[j]["m_jj"+drpt], "m_jj"+dr+drpt, 80, 0, 4000);
	  }
	}
	
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      for (size_t j = 0; j < DRS.size(); ++j) {
	std::string dr = "_xxx";
	sprintf(&dr[1], "%0.1f", DRS[j]);
	// Get Higgs
	const Particles higgses = apply<ParticleFinder>(event, "Higgses").particles();
	if (higgses.empty()) continue;
	if (higgses.size() > 1) continue;
	const Particle higgs = higgses[0];
	const double ptH = higgs.pT();
	const double yh  = higgs.rap();

	// Get jets and leading dijet system
	const Jets jets = apply<FastJets>(event, "Jets"+dr).jetsByPt(Cuts::pT > 30*GeV && Cuts::absrap < 4.4);
	if (jets.size() < 2) continue;
	const int njets = jets.size();
	const double htj = sum(jets, Kin::pT, 0.0);
	const double ht = htj + higgs.pT();

	// Get leading dijet system
	const FourMomentum dijet = jets[0].mom() + jets[1].mom();
	const double m12 = dijet.mass();
	const double pt12 = dijet.pT();
	const double y12 = dijet.rap();
	const double dy12 = std::abs(jets[0].rap() - jets[1].rap());
	const double dphi12 = deltaPhi(jets[0], jets[1]);
	const double dR12   = sqrt(dy12*dy12+dphi12*dphi12);
	const double y1  = jets[0].rap();
	const double y2  = jets[1].rap();

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
	  _c_xs[j][i]->fill();
	  _h[j][i]["deltaphi_jj"]->fill(dphi12);
	  _h[j][i]["y_j12_02bin"]->fill(y12);
	  _h[j][i]["pth_largebin"]->fill(ptH/GeV);
	  _h[j][i]["pth_finebin"]->fill(ptH/GeV);
	  _h[j][i]["njets"]->fill(njets);
	  _h[j][i]["HT"]->fill(ht/GeV);
	  _h[j][i]["pthjj_largebin"]->fill(ptjjH/GeV);
	  _h[j][i]["logpthjj_largebin"]->fill(log10(ptjjH/GeV));
	  _h[j][i]["ptjj_largebin"]->fill(pt12/GeV);
	  /*for (size_t c = 0; c < HPTCUTS.size(); ++c) {
	    if (ptH/GeV>HPTCUTS[c]) {
	    std::string pt = "_xxx";
	    sprintf(&pt[1], "%0.0f", HPTCUTS[c]);
	    _h[j][i]["deltay_jj"+pt]->fill(dy12);
	    _h[j][i]["deltaR_jj"+pt]->fill(sqrt(dy12*dy12+dphi12*dphi12));
	    _h[j][i]["m_jj"+pt]->fill(m12);
	    }
	    }*/
	}

	// Fill dijet mass histograms
	_h_mjj[j]["mjj_STXS"]->fill(m12/GeV);
	_h_mjj[j]["mjj_ATLAS_CONF_029"]->fill(m12/GeV);
	_h_mjj[j]["mjj_100GeVbin"]->fill(m12/GeV);
	if (dy12 < 2.)
	  _h_mjj[j]["mjj_100GeVbin_center"]->fill(m12/GeV);
	if (dy12 >= 2. && dy12 < 4.)
	  _h_mjj[j]["mjj_100GeVbin_middle"]->fill(m12/GeV);
	if (dy12 > 4.)
	  _h_mjj[j]["mjj_100GeVbin_forward"]->fill(m12/GeV);
	if (dy12 > 3 && m12 > 400*GeV && dphijjH > 2.8)
	  _h_mjj[j]["mjj_100GeVbin_ATLAS"]->fill(m12/GeV);

	// // Now the dy12 histograms, with various m12 cuts
	// _h_dyjj[j]["deltay_jj"]->fill(dy12);
	// if (m12 < 350*GeV)
	//   _h_dyjj[j]["deltay_jj_light"]->fill(dy12);
	// else // if (m12 < 350*GeV)
	//   _h_dyjj[j]["deltay_jj_heavy"]->fill(dy12);

	for (size_t c = 0; c < HPTCUTS.size(); ++c) {
	  if (ptH/GeV>HPTCUTS[c]) {
	    std::string pt = std::to_string(int(HPTCUTS[c]));
	    pt = "_pt"+pt;
	    _hh[j]["deltay_jj"+pt]->fill(dy12);
	    _hh[j]["deltaR_jj"+pt]->fill(dR12);
	    _hh[j]["deltaphi_jj"+pt]->fill(dphi12);
	    _hh[j]["y_h"+pt]->fill(yh);
	    _hh[j]["y_j1"+pt]->fill(y1);
	    _hh[j]["y_j2"+pt]->fill(y2);
	    if(jets.size()>2) _hh[j]["y_j3"+pt]->fill(jets[2].rap());
	    _hh[j]["m_jj"+pt]->fill(m12);

	    
	  }
	}
	
	for (size_t c = 0; c < HPTCUTSV2.size(); ++c) {
	  if (ptH/GeV>HPTCUTSV2[c]) {
	    std::string pt = std::to_string(int(HPTCUTSV2[c]));
	    pt = "_pt"+pt;
	    for (size_t c1 = 0; c1 < M12MAXCUTS.size(); ++c1) {
	      if(m12/GeV < M12MAXCUTS[c1]){
		std::string mpt = std::to_string(int(M12MAXCUTS[c1]));
		mpt = "_m"+mpt+pt;
		_hh[j]["deltay_jj"+mpt]->fill(dy12);
		_hh[j]["deltaR_jj"+mpt]->fill(dR12);
		_hh[j]["deltaphi_jj"+mpt]->fill(dphi12);
		_hh[j]["y_h"+mpt]->fill(yh);
		_hh[j]["y_j1"+mpt]->fill(y1);
		_hh[j]["y_j2"+mpt]->fill(y2);
		if(jets.size()>2) _hh[j]["y_j3"+mpt]->fill(jets[2].rap());
	      }
	    }
	  }
	}

	for (size_t c = 0; c < HPTCUTSV2.size(); ++c) {
	  if (ptH/GeV>HPTCUTSV2[c]) {
	    std::string pt = std::to_string(int(HPTCUTSV2[c]));
	    pt = "_pt"+pt;
	    for (size_t c1 = 0; c1 < DR12CUTS.size(); ++c1) {
	      if(dR12 < DR12CUTS[c1]){
		std::string drpt =  "_drxxx";
		sprintf(&drpt[3], "%0.1f", DR12CUTS[c1]);
		drpt = drpt+pt;
		_hh[j]["m_jj"+drpt]->fill(m12);
	      }
	    }
	  }
	}
	
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for (size_t j = 0; j < DRS.size(); ++j) {
	for (size_t i = 0; i < SELNAMES.size(); ++i) {
	  scale(_c_xs[j][i], crossSection()/femtobarn/sumOfWeights());
	  for (auto kv : _h[j][i])
	    scale(kv.second, crossSection()/femtobarn/sumOfWeights());
	}
	for (auto kv : _h_mjj[j])
	  scale(kv.second, crossSection()/femtobarn/sumOfWeights());
	for (auto kv : _hh[j])
	  scale(kv.second, crossSection()/femtobarn/sumOfWeights());
	
	/* for (auto kv : _h_dyjj[j])
	   scale(kv.second, crossSection()/femtobarn/sumOfWeights());*/
      }
    }


    /// @name Histograms
    //@{
    map<string, Histo1DPtr> _h[15][8], _hh[15];
    CounterPtr _c_xs[15][8];
    map<string,Histo1DPtr> _h_mjj[15]; //, _h_dyjj[15];
    //@}

    static const vector<double> DRS, HPTCUTS;

    /// Cut values for standard histogram sets (other than m12 and dy12)
    //@{
    static const vector<string> SELNAMES;
    static const vector<doubles> M12CUTS;
    static const vector<doubles> DY12CUTS;
    static const vector<double> DPHIJHCUTS;
    static const vector<double> HPTCUTSV2, M12MAXCUTS, DR12CUTS;
    //@}

  };


  // Static const initializers
  //const vector<double>  MC_HJETSVBF::DRS = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5 };
  const vector<double>  MC_HJETSVBF::DRS = {0.2, 0.4, 0.7, 1.0};
  const vector<string>  MC_HJETSVBF::SELNAMES = {"", "_light_center", "_heavy_center", "_light_middle", "_heavy_middle", "_light_forward", "_heavy_forward", "_ATLAS"};
  const vector<doubles> MC_HJETSVBF::M12CUTS  = {{0., HUGE_VAL}, {0., 350.}, {350., HUGE_VAL}, {0., 350.}, {350., HUGE_VAL}, {0., 350.}, {350., HUGE_VAL}, {400., HUGE_VAL}};
  const vector<doubles> MC_HJETSVBF::DY12CUTS = {{0., HUGE_VAL}, {0., 2.}, {0., 2.}, {2., 4.}, {2., 4.}, {4., HUGE_VAL}, {4., HUGE_VAL}, {3., HUGE_VAL}};
  const vector<double>  MC_HJETSVBF::DPHIJHCUTS = {0., 0., 0., 0., 0., 0., 0., 2.8};
  const vector<double> MC_HJETSVBF::HPTCUTS  = {0., 100.,200.,300.,400.,500.,1000.};
  const vector<double> MC_HJETSVBF::HPTCUTSV2= {0., 200., 500.};
  const vector<double> MC_HJETSVBF::M12MAXCUTS={100.,200.,300.,400.,500.,1000.,2000.};
  const vector<double> MC_HJETSVBF::DR12CUTS={0.5, 1.0, 1.5,  2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0};
  DECLARE_RIVET_PLUGIN(MC_HJETSVBF);

}
