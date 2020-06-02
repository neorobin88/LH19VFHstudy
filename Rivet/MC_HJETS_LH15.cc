// -*- C++ -*-

// ChangeLog
// -----------------------------------------------------------------------------
// 02/06
// * filling underflow and overflow bins
// -----------------------------------------------------------------------------
// 13/06
// * bugfix
// -----------------------------------------------------------------------------
// 15/05
// * added some WBF and WBF2 observables
// -----------------------------------------------------------------------------
// 19/05 MS
// * changed name to MC_HJETS_LH15
// * switched all histograms to NLO histograms to properly do fixed-order NLO
//   statistical uncertainties
// * make analysis work with h->yy decays and stable higgs final state alike
// * remove all fiducial cuts
// * fix RapJets definition
// * fix sum_tau_jet histogramming
// * added exclusive jet multiplicities where missing
// * commented needless filling of underflow and overflow bins
// * change deltaphi_jj_bins to be in multiples of pi
// * added 3j observables (pT(H),pT(j) incl. and excl.)
// * grouped histogram initialisation to be readable/maintainable
// * used multiple bins to book-keep loose and tight cross sections
//   (0..dijet, 1..VBF, 2..VBF2)
// * some code cosmetics
// -----------------------------------------------------------------------------
// 21/05
// * bugfix in tight cross section histogramming
// * bugfix in histogram synchronisation
// -----------------------------------------------------------------------------
// 03/06
// * renamed some observables sensibly
// -----------------------------------------------------------------------------
// 05/06
// * add dR(y,j1) and dR(y,j2)
// * add jet veto cross sections
// -----------------------------------------------------------------------------
// 06/06
// * remove coarsly binned histograms
// * adjusted binnings of many others to half-way match statistical population
// -----------------------------------------------------------------------------
// 09/07
// * make FastJets object a member to fix memleak
// -----------------------------------------------------------------------------
// 28/10
// * fix jet1_y, jet3_y
// -----------------------------------------------------------------------------
// 19/11
// * bugfix mindy
// -----------------------------------------------------------------------------
// 01/06/2020 - Andy Buckley
// * convert to Rivet 3


#include "Rivet/Analyses/MC_JetAnalysis.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"

namespace Rivet {


  class MC_HJETS_LH15 : public Analysis {
  private:

    double _jrap, _jR, _jpT;
    double _mH, _mHdev;
    double _wbfdyjj, _wbfmjj;
    double _dphiHjjtight;
    double _taujcut;
    std::map<std::string,Histo1DPtr> histos;

  public:


    MC_HJETS_LH15() :
      Analysis("MC_HJETS_LH15"),
      _jrap(4.4), _jR(0.4), _jpT(30.*GeV), _mH(125.*GeV), _mHdev(1.*GeV),
      _wbfdyjj(2.8), _wbfmjj(400.*GeV), _dphiHjjtight(2.6), _taujcut(8.*GeV)
    {}

    // double sqr(const double& x) { return x*x; }

    // Will be used with p1=(1,0,0,1),p3=(1,0,0,-1). Only non-zero terms kept.
    std::complex<double> EPSTENSOR(const FourMomentum &p1,
                                   const FourMomentum &p2,
                                   const FourMomentum &p3,
                                   const FourMomentum &p4)
    {
      return -std::complex<double>(0.,1.)
        *(-p1.z()*p2.x()*p3.E()*p4.y()+p1.z()*p2.y()*p3.E()*p4.x()
          +p1.E()*p2.x()*p3.z()*p4.y()-p1.E()*p2.y()*p3.z()*p4.x());
    }

    void fillVetoCrossSection(double pT, const string& id) {
      Histo1DPtr jvh = histos[id];
      size_t index = jvh->binIndexAt(pT);
      for (size_t i = index; i < jvh->numBins(); ++i) {
        jvh->fillBin(i, jvh->bin(i).xWidth());
      }
    }

    void init() {
      FinalState fs;
      IdentifiedFinalState higgses(PID::HIGGS);
      IdentifiedFinalState photons(PID::PHOTON);
      VetoedFinalState rest(fs);
      rest.addVetoOnThisFinalState(higgses);
      rest.addVetoOnThisFinalState(photons);
      declare(fs, "FS");
      declare(higgses, "Higgses");
      declare(photons, "Photons");
      declare(rest, "Rest");
      declare(FastJets(rest, FastJets::ANTIKT, _jR), "Jets");

      inithistos();
    }


    void inithistos() {

      book(histos["XS"], "XS",1,0.,1.);
      book(histos["m_gammagamma"], "m_gammagamma", bwspace(21,124.99,125.01,125.,0.00407));

      const vector<double> H_pT_bins{0.,5.,10.,15.,20.,25.,30.,35.,40.,45.,50.,55.,60.,
          65.,70.,75.,80.,85.,90.,95.,100.,110.,120.,130.,140.,150.,
          160.,170.,180.,190.,200.};

      const vector<double> H_pT_jj_bins{0.,10.,20.,30.,40.,50.,60.,70.,80.,90.,100.,110.,120.,
          130.,140.,150.,160.,170.,180.,190.,200.,220.,240.,260.,
          280.,300.,320.,340.,360.,380.,400.};

      const vector<double> Hj_pT_bins{0.,5.,10.,15.,20.,25.,30.,35.,40.,45.,50.,55.,60.,
          65.,70.,75.,80.,85.,90.,95.,100.,110.,120.,130.,140.,150.,
          160.,170.,180.,190.,200.,220.,240.,260.,280.,300.,
          320.,340.,360.,380.,400.};

      const vector<double> H_pT_0j_excl_bins{0.,20.,30.,45.,200.};
      const vector<double> H_pT_1j_excl_bins{0.,40.,60.,95.,200.};
      const vector<double> H_pT_2j_excl_bins{0.,90.,140.,200.};
      const vector<double> H_pT_3j_excl_bins{0.,90.,140.,200.};
      const vector<double> H_pT_0j_incl_bins{0.,20.,30.,45.,200.};
      const vector<double> H_pT_1j_incl_bins{0.,40.,60.,95.,200.};
      const vector<double> H_pT_2j_incl_bins{0.,90.,140.,200.};
      const vector<double> H_pT_3j_incl_bins{0.,90.,140.,200.};

      const vector<double> jet1_pT_bins{0.,30.,50.,70.,100.,140.,500.};
      const vector<double> jet2_pT_bins{0,30,40,50,140,500.};
      const vector<double> jet3_pT_bins{0,30,50,150,500.};

      const vector<double> H_y_bins{0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,
          1.4,1.5,1.6,1.7,1.8,1.9,2.,2.1,2.2,2.3,2.4,2.5,2.7,2.9,
          3.1,3.3,3.5,4.0,4.5,5};

      const vector<double> jet_y_bins{0.,0.25,0.5,0.75,1.,1.25,1.5,1.75,2.,2.25,2.5,
          3.,3.5,4.,4.4};

      const vector<double> deltaphi_jj_bins{0.,PI/16.,2.*PI/16.,3.*PI/16.,4.*PI/16.,
          5.*PI/16.,6.*PI/16.,7.*PI/16.,8.*PI/16.,
          9.*PI/16.,10.*PI/16.,11.*PI/16.,12.*PI/16.,
          13.*PI/16.,14.*PI/16.,15.*PI/16.,PI};

      const vector<double> deltay_jj_bins{0.,0.5,1.,1.5,2.,2.5,3.,3.5,4.,4.5,5.,5.5,6.,
          6.5,7.,7.5,8.,8.5,9.,9.5,10.};

      const vector<double> dijet_mass_bins{0.,50.,100.,150.,200.,250.,300.,350.,400.,450.,
          500.,550.,600.,650.,700.,750.,800.,850.,900.,
          950.,1000.};
      const vector<double> deltay_yy_bins{0,0.3,0.6,0.9,1.2,1.5,2.0,2.55};
      const vector<double> dR_y_j_bins{0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.2,1.4,1.6,1.8,
          2.0,2.4,2.8,3.2,3.6,4.0,4.5,5.0,5.5,6.0,7.0,8.0,10.0};

      const vector<double> tau_jet_bins{0.,4.,8.,12.,16.,20.,25.,30.,40.,60.,85.};
      const vector<double> HT_bins{0.,30.,40.,50.,60.,70.,90.,110.,130.,150.,200.,
          250.,300.,400.,500.,600.,800.,1000.};
      const vector<double> pTt_bins{0.,10.,20.,30.,40.,60.,80.,150.,500.};
      const vector<double> deltaphi_Hjj_bins{0.0,1.0,2.0,2.3,2.6,2.8,2.9,3.,3.05,3.1,PI};

      const vector<double> deltay_H_jj_bins{0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,
          2.6,2.8,3.0,3.5,4.0,5.0,6.0,8.0};

      const vector<double> delta_phi2_bins{-PI,-3.0,-2.8,-2.6,-2.4,-2.2,-2.0,-1.8,-1.6,-1.4,-1.2,
          -1.,-0.8,-0.6,-0.4,-0.2,0.,0.2,0.4,0.6,0.8,1.,1.2,1.4,
          1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,PI};

      const vector<double> H_dijet_mass_bins{0.,100.,150.,200.,250.,300.,350.,400.,450.,500.,550.,
          600.,650.,700.,750.,800.,850.,900.,950.,1000.,1100.,
          1200.,1300.,1400.,1500.,1700.,2000.};

      // fine binnings always end on ""
      // incl. and excl. jet multis
      book(histos["NJet_incl_30"], "NJet_incl_30",4,-0.5,3.5);
      book(histos["NJet_excl_30"], "NJet_excl_30",4,-0.5,3.5);
      book(histos["NJet_incl_50"], "NJet_incl_50",4,-0.5,3.5);
      book(histos["NJet_excl_50"], "NJet_excl_50",4,-0.5,3.5);
      book(histos["NJet_incl_30_jhj"], "NJet_incl_30_jhj",4,-0.5,3.5);
      book(histos["NJet_excl_30_jhj"], "NJet_excl_30_jhj",4,-0.5,3.5);
      book(histos["NJet_excl_30_VBF"], "NJet_excl_30_VBF",4,-0.5,3.5);
      book(histos["NJet_incl_30_VBF"], "NJet_incl_30_VBF",4,-0.5,3.5);
      book(histos["NJet_excl_30_VBF2"], "NJet_excl_30_VBF2",4,-0.5,3.5);
      book(histos["NJet_incl_30_VBF2"], "NJet_incl_30_VBF2",4,-0.5,3.5);

      // pT(H) in incl. and excl. jet bins
      book(histos["H_pT_incl"], "H_pT_incl",H_pT_bins);
      book(histos["H_pT_excl"], "H_pT_excl",H_pT_bins);

      book(histos["H_j_pT_incl"], "H_j_pT_incl",H_pT_bins);
      book(histos["H_j_pT_excl"], "H_j_pT_excl",H_pT_bins);

      book(histos["H_jj_pT_incl"], "H_jj_pT_incl",H_pT_jj_bins);
      book(histos["H_jj_pT_excl"], "H_jj_pT_excl",H_pT_jj_bins);

      book(histos["H_jjj_pT_incl"], "H_jjj_pT_incl",H_pT_jj_bins);
      book(histos["H_jjj_pT_excl"], "H_jjj_pT_excl",H_pT_jj_bins);

      book(histos["H_jj_pT_VBF"], "H_jj_pT_VBF",H_pT_jj_bins);
      book(histos["H_jj_pT_VBF2"], "H_jj_pT_VBF2",H_pT_jj_bins);

      // pT(H+nj) in incl. and excl. jet bins
      book(histos["Hj_pT_incl"], "Hj_pT_incl",Hj_pT_bins);
      book(histos["Hj_pT_excl"], "Hj_pT_excl",Hj_pT_bins);

      book(histos["Hjj_pT_incl"], "Hjj_pT_incl",Hj_pT_bins);
      book(histos["Hjj_pT_excl"], "Hjj_pT_excl",Hj_pT_bins);

      // pT(j) in incl. and excl. jet bins
      book(histos["jet1_pT_incl"], "jet1_pT_incl",H_pT_bins);
      book(histos["jet1_pT_excl"], "jet1_pT_excl",H_pT_bins);

      book(histos["jet2_pT_incl"], "jet2_pT_incl",H_pT_bins);
      book(histos["jet2_pT_excl"], "jet2_pT_excl",H_pT_bins);

      book(histos["jet3_pT_incl"], "jet3_pT_incl",H_pT_bins);
      book(histos["jet3_pT_excl"], "jet3_pT_excl",H_pT_bins);

      // inclusive Higgs and jet rapidities
      book(histos["H_y"], "H_y",H_y_bins);

      book(histos["jet1_y"], "jet1_y",jet_y_bins);

      book(histos["jet2_y"], "jet2_y",jet_y_bins);

      book(histos["jet3_y"], "jet3_y",jet_y_bins);

      // photon observables
      book(histos["cos_theta_star"], "cos_theta_star",10,0.,1.);
      book(histos["cos_theta_star_80"], "cos_theta_star_80",4,0.,1.);
      book(histos["cos_theta_star_200"], "cos_theta_star_200",4,0.,1.);
      book(histos["cos_theta_star_gt200"], "cos_theta_star_gt200",4,0.,1.);
      book(histos["deltay_yy"], "deltay_yy",deltay_yy_bins);
      book(histos["dR_y_j1"], "dR_y_j1",dR_y_j_bins);
      book(histos["dR_y_j2"], "dR_y_j2",dR_y_j_bins);

      // book-keep loose and tight cross sections
      book(histos["loose"], "loose",3,-0.5,2.5);
      book(histos["tight"], "tight",3,-0.5,2.5);

      // \Delta\phi(jj) in incl. and excl. jet bins
      book(histos["deltaphi_jj_incl"], "deltaphi_jj_incl",deltaphi_jj_bins);
      book(histos["deltaphi_jj_excl"], "deltaphi_jj_excl",deltaphi_jj_bins);
      book(histos["deltaphi_jj_VBF"], "deltaphi_jj_VBF",deltaphi_jj_bins);
      book(histos["deltaphi_jj_VBF2"], "deltaphi_jj_VBF2",deltaphi_jj_bins);

      // \Delta y(jj)
      book(histos["deltay_jj"], "deltay_jj",deltay_jj_bins);

      // m(jj)
      book(histos["dijet_mass"], "dijet_mass",dijet_mass_bins);

      // m(Hjj)
      book(histos["H_dijet_mass"], "H_dijet_mass",H_dijet_mass_bins);

      // \Delta\phi(H,jj) incl. and excl.
      book(histos["deltaphi_Hjj_incl"], "deltaphi_Hjj_incl",deltaphi_Hjj_bins);
      book(histos["deltaphi_Hjj_excl"], "deltaphi_Hjj_excl",deltaphi_Hjj_bins);
      book(histos["deltaphi_Hjj_VBF"], "deltaphi_Hjj_VBF",deltaphi_Hjj_bins);
      book(histos["deltaphi_Hjj_VBF2"], "deltaphi_Hjj_VBF2",deltaphi_Hjj_bins);

      // \Delta y(H,jj)
      book(histos["deltay_H_jj"], "deltay_H_jj",deltay_H_jj_bins);

      // HT
      book(histos["HT_all"], "HT_all",HT_bins);
      book(histos["HT_jets"], "HT_jets",HT_bins);

      // tau(j) observables
      book(histos["tau_jet1"], "tau_jet1",tau_jet_bins);
      book(histos["tau_jet2"], "tau_jet2",tau_jet_bins);
      book(histos["tau_jet3"], "tau_jet3",tau_jet_bins);
      book(histos["tau_jet_max"] , "tau_jet_max",tau_jet_bins);
      book(histos["sum_tau_jet"] , "sum_tau_jet",tau_jet_bins);

      // pTt
      book(histos["pTt"], "pTt",pTt_bins);

      // \phi_2
      book(histos["deltaphi2"], "deltaphi2",delta_phi2_bins);
      book(histos["deltaphi2_VBF"], "deltaphi2_VBF",delta_phi2_bins);
      book(histos["deltaphi2_VBF2"], "deltaphi2_VBF2",delta_phi2_bins);


      // IVAN *****************************************************

      string jj[] = { "pT", "dy" };

      for (int i=0; i<2; ++i) {
        book(histos["jj"+jj[i]+"_dy"], "jj"+jj[i]+"_dy",18,0,9);
        book(histos["jj"+jj[i]+"_dy_2j_excl"], "jj"+jj[i]+"_dy_2j_excl",18,0,9);
        book(histos["jj"+jj[i]+"_dy_3j_excl"], "jj"+jj[i]+"_dy_3j_excl",18,0,9);
      }

      // rapidity distance between forward and backward jets
      book(histos["jjfb_dy"], "jjfb_dy",18,0,9);

      for (int i=0; i<2; ++i) {
        for (int j=1; j<=3; ++j) {
          for (int y=1; y<=6; ++y) {
            stringstream ss;
            ss << "jet" << j << "_pT_jj" << jj[i] << "_mindy" << y;
            book(histos[ss.str()], ss.str(),30,0,300);
          }
        }
      }

      book(histos["xs_central_jet_veto"], "xs_central_jet_veto",25,0.,5.);
      book(histos["xs_central_jet_veto_VBF"], "xs_central_jet_veto_VBF",25,0.,5.);
      book(histos["xs_central_jet_veto_VBF2"], "xs_central_jet_veto_VBF2",25,0.,5.);


      // END IVAN *************************************************
      book(histos["xs_jet_veto_j0"], "xs_jet_veto_j0",logspace(300,1.,1000.));
      book(histos["xs_jet_veto_j1_30"], "xs_jet_veto_j1_30",logspace(300,1.,1000.));
      book(histos["xs_jet_veto_j1_50"], "xs_jet_veto_j1_50",logspace(300,1.,1000.));
      book(histos["xs_jet_veto_j1_100"], "xs_jet_veto_j1_100",logspace(300,1.,1000.));
      book(histos["xs_jet_veto_j1_200"], "xs_jet_veto_j1_200",logspace(300,1.,1000.));
      book(histos["xs_jet_veto_j1_500"], "xs_jet_veto_j1_500",logspace(300,1.,1000.));
      book(histos["xs_jet_veto_j1_1000"], "xs_jet_veto_j1_1000",logspace(300,1.,1000.));
      book(histos["xs_jet_veto_h_50"], "xs_jet_veto_h_50",logspace(300,1.,1000.));
      book(histos["xs_jet_veto_h_100"], "xs_jet_veto_h_100",logspace(300,1.,1000.));
      book(histos["xs_jet_veto_h_200"], "xs_jet_veto_h_200",logspace(300,1.,1000.));
      book(histos["xs_jet_veto_h_500"], "xs_jet_veto_h_500",logspace(300,1.,1000.));
    }


    /// Do the analysis
    void analyze(const Event & e) {

      Particles higgses = apply<IdentifiedFinalState>(e, "Higgses").particles();
      Particles photons = apply<IdentifiedFinalState>(e, "Photons").particles();
      Particles rest = apply<VetoedFinalState>(e, "Rest").particles();

      // require either one stable Higgs or at least two photons
      if (higgses.size()>1) vetoEvent;
      FourMomentum hmom;
      size_t idph1(0),idph2(0);
      vector<FourMomentum> phs;
      if (higgses.size()==1) {
        hmom = higgses[0].momentum();
      }
      else if (photons.size()>1) {
        // reconstruct Higgs from photon pair with correct inv. mass
        // only take first one
        bool foundone(false);
        for (size_t i(0);i<photons.size();++i) {
          for (size_t j(i+1);j<photons.size();++j) {
            if (!foundone &&
                fabs((photons[i].momentum()
                      +photons[j].momentum()).mass()-_mH)<_mHdev) {
              idph1=i; idph2=j;
              hmom = photons[i].momentum()+photons[j].momentum();
              phs.push_back(photons[i].momentum());
              phs.push_back(photons[j].momentum());
              foundone=true;
              break;
            }
          }
          if (foundone) break;
        }
      }
      else vetoEvent;

      // check that found one Higgs
      if (higgses.size()==0 && phs.size()!=2) vetoEvent;

      // add remaining photons to the remaining final state
      for (size_t i(0); i<photons.size(); ++i) {
        if (idph1==idph2 || (i!=idph1 && i!=idph2)) rest.push_back(photons[i]);
      }

      // calculate jets
      //_jetalgo.calc(rest);

      // Create y ordered and pt ordered jet vectors
      Jets PTJets,RapJets,alljets;

      const PseudoJets jets = apply<FastJets>(e, "Jets").pseudoJetsByPt(0.*GeV);
      for (const Jet& jetcand : jets) {
        if (fabs(jetcand.momentum().rapidity()) < _jrap) {
          alljets.push_back(jetcand);
        }
      }
      for (const Jet& jetcand : jets) {
        if (fabs(jetcand.momentum().rapidity()) < _jrap) {
          PTJets.push_back(jetcand);
        }
      }
      for (const Jet& jetcand : jets) {
        if (fabs(jetcand.momentum().rapidity()) < _jrap) {
          RapJets.push_back(jetcand);
        }
      }

      // cross check
      if (PTJets.size()!=RapJets.size()) abort();

      histos["XS"]->fill(0.5);

      // photon observables
      if (phs.size()>1) {
        histos["m_gammagamma"]->fill(hmom.mass()/GeV);
        // |costheta*| from 1307.1432
        double cts(abs(sinh(phs[0].eta()-phs[1].eta()))/
                   sqrt(1.+sqr(hmom.pT()/hmom.mass()))
                   * 2.*phs[0].pT()*phs[1].pT()/sqr(hmom.mass()));
        histos["cos_theta_star"]->fill(cts);

        if (hmom.pT()<80.*GeV) {
          histos["cos_theta_star_80"]->fill(cts);
        }

        if (hmom.pT()>80.*GeV && hmom.pT()<200.*GeV) {
          histos["cos_theta_star_200"]->fill(cts);
        }

        if (hmom.pT()>200.*GeV) {
          histos["cos_theta_star_gt200"]->fill(cts);
        }

        double pTt = fabs(phs[0].px()*phs[1].py()-phs[1].px()*phs[0].py())/
          ((phs[0]-phs[1]).pT()*2);
        histos["pTt"]->fill(pTt);
        double deltay_yy = fabs(phs[0].rapidity()-phs[1].rapidity());
        histos["deltay_yy"]->fill(deltay_yy);
      }

      // inclusive histograms

      histos["H_pT_incl"]->fill(hmom.pT()/GeV);

      histos["H_y"]->fill(fabs(hmom.rapidity()));

      histos["NJet_excl_30"]->fill(PTJets.size());
      for (size_t i(0);i<4;++i) {
        if (PTJets.size()>=i) histos["NJet_incl_30"]->fill(i);
      }
      histos["NJet_incl_50"]->fill(0);

      // 0j jet veto cross section
      double jv0pT1 = alljets.size()>0?alljets[0].momentum().pT():0.;
      fillVetoCrossSection(jv0pT1,"xs_jet_veto_j0");

      // njets == 0;
      if (PTJets.size()==0) {
        histos["H_pT_excl"]->fill(hmom.pT()/GeV);
        // why????
        // 6/2 added fill for jet1_pT for 0-30 GeV bin, i.e. no jets
        // histos["jet1_pT_incl"]->fill(10);
        // 6/2 added fill for overflow bin for 0 jets in event
        // histos["deltay_jj"]->fill(9);
        // 6/2 added fill for overflow bin for 0 jets in event
        // histos["Hjj_pT_incl"]->fill(160);
        // 6/2 added fill for overflow bin for 0 jets in event
        // histos["jet2_y"]->fill(4.6);
        // 6/2 added fill for overflow bin for 0 jets in event
        // histos["jet2_pT_incl"]->fill(400);
      }

      // njets > 0;
      if (PTJets.size()>0) {
        const FourMomentum& j1(PTJets[0].momentum());

        histos["jet1_pT_incl"]->fill(j1.pT()/GeV);
        histos["jet1_y"]->fill(fabs(j1.rapidity()));
        histos["Hj_pT_incl"]->fill((hmom+j1).pT()/GeV);
        histos["H_j_pT_incl"]->fill(hmom.pT()/GeV);

        // Calculate tau
        double tauJet1 = sqrt(sqr(j1.pT()) + sqr(j1.mass()))/
          (2.*cosh(j1.rapidity() - hmom.rapidity()));

        histos["tau_jet1"]->fill(tauJet1/GeV);

        // njets == 1;
        if (PTJets.size()==1) {
          histos["Hj_pT_excl"]->fill((hmom+j1).pT()/GeV);
          histos["H_j_pT_excl"]->fill(hmom.pT()/GeV);
          histos["jet1_pT_excl"]->fill(j1.pT()/GeV);
          // again, why???
          // 6/2 added fill for j2_pT for 0-30 GeV bins, i.e. no 2nd jet
          // histos["jet2_pT_incl"]->fill(10);
          // 6/2 added fill for overflow bin for 1 jet in event
          // histos["deltay_jj"]->fill(9);
          // 6/2 added fill for overflow bin for 1 jet in event
          // histos["Hjj_pT_incl"]->fill(160);
          // 6/2 added fill for overflow bin for 1 jet in event
          // histos["jet2_y"]->fill(4.6);
        }

        // 1j jet veto cross section with minimal pT(j1)
        double jv1pT2 = alljets.size()>1?alljets[1].momentum().pT():0.;
        if (alljets[0].momentum().pT()>30.*GeV)
          fillVetoCrossSection(jv1pT2,"xs_jet_veto_j1_30");
        if (alljets[0].momentum().pT()>50.*GeV)
          fillVetoCrossSection(jv1pT2,"xs_jet_veto_j1_50");
        if (alljets[0].momentum().pT()>100.*GeV)
          fillVetoCrossSection(jv1pT2,"xs_jet_veto_j1_100");
        if (alljets[0].momentum().pT()>200.*GeV)
          fillVetoCrossSection(jv1pT2,"xs_jet_veto_j1_200");
        if (alljets[0].momentum().pT()>500.*GeV)
          fillVetoCrossSection(jv1pT2,"xs_jet_veto_j1_500");
        if (alljets[0].momentum().pT()>1000.*GeV)
          fillVetoCrossSection(jv1pT2,"xs_jet_veto_j1_1000");
      }

      // 1j jet veto cross section with minimal pT(j1)
      double jv1pT2 = alljets.size()>1?alljets[1].momentum().pT():0.;
      if (hmom.pT()>50.*GeV)
        fillVetoCrossSection(jv1pT2,"xs_jet_veto_h_50");
      if (hmom.pT()>100.*GeV)
        fillVetoCrossSection(jv1pT2,"xs_jet_veto_h_100");
      if (hmom.pT()>200.*GeV)
        fillVetoCrossSection(jv1pT2,"xs_jet_veto_h_200");
      if (hmom.pT()>500.*GeV)
        fillVetoCrossSection(jv1pT2,"xs_jet_veto_h_500");

      // njets > 1;
      if (PTJets.size()>1) {
        const FourMomentum& j1(PTJets[0].momentum());
        const FourMomentum& j2(PTJets[1].momentum());

        if (phs.size()>1) {
          // fill dR of j1 and j2 and nearest photon
          histos["dR_y_j1"]->fill(min(deltaR(j1,phs[0]),deltaR(j1,phs[1])));
          histos["dR_y_j2"]->fill(min(deltaR(j2,phs[0]),deltaR(j2,phs[1])));
        }

        const FourMomentum& jb(RapJets.front().momentum());
        const FourMomentum& jf(RapJets.back().momentum());

        // Calculation of phi_2 from arXiv:1001.3822
        // phi_2 = azimuthal angle between the vector sum of jets
        // forward and jets backward of the Higgs boson

        FourMomentum vsumf(0.,0.,0.,0.);
        FourMomentum vsumb(0.,0.,0.,0.);
        FourMomentum p1(1.,0.,0.,1.);
        FourMomentum p2(1.,0.,0.,-1.);

        bool f_nonzero(false),b_nonzero(false);
        for (const Jet& jj : RapJets) {
          if (jj.momentum().rapidity()>hmom.rapidity()) {
            vsumf += jj.momentum();
            f_nonzero = true;
          }
          else {
            vsumb += jj.momentum();
            b_nonzero = true;
          }
        }

        double phi2(-10.);
        // Calculate phi_2
        if (f_nonzero && b_nonzero) {
          phi2 = acos((vsumb.x()*vsumf.x()+vsumb.y()*vsumf.y())/
                      (sqrt(vsumb.x()*vsumb.x()+vsumb.y()*vsumb.y())*
                       sqrt(vsumf.x()*vsumf.x()+vsumf.y()*vsumf.y())));
          if (imag(EPSTENSOR(p1,vsumb,p2,vsumf))<0.) phi2 *= -1.;

        }
        else if (!f_nonzero) {
          vsumb -= jf;
          phi2 = acos((vsumb.x()*jf.x()+vsumb.y()*jf.y())/
                      (sqrt(vsumb.x()*vsumb.x()+vsumb.y()*vsumb.y())*
                       sqrt(jf.x()*jf.x()+jf.y()*jf.y())));
          if (imag(EPSTENSOR(p1,vsumb,p2,jf))<0.)    phi2 *= -1.;
        }
        else {
          vsumf -= jb;
          phi2 = acos((jb.x()*vsumf.x()+jb.y()*vsumf.y())/
                      (sqrt(jb.x()*jb.x()+jb.y()*jb.y())*
                       sqrt(vsumf.x()*vsumf.x()+vsumf.y()*vsumf.y())));
          if (imag(EPSTENSOR(p1,jb,p2,vsumf))<0.)    phi2 *= -1.;
        }
        histos["deltaphi2"]->fill(phi2);

        if (f_nonzero && b_nonzero) {
          histos["NJet_excl_30_jhj"]->fill(PTJets.size());
          for (size_t i(0);i<4;++i) {
            if (PTJets.size()>=i) histos["NJet_incl_30_jhj"]->fill(i);
          }
        }

        histos["deltaphi_jj_incl"]->fill(deltaPhi(j1,j2));
        histos["deltaphi_Hjj_incl"]->fill(deltaPhi(hmom,j1+j2));
        histos["Hjj_pT_incl"]->fill((hmom+j1+j2).pT()/GeV);
        histos["H_jj_pT_incl"]->fill(hmom.pT()/GeV);
        histos["jet2_pT_incl"]->fill(j2.pT()/GeV);
        histos["jet2_y"]->fill(fabs(j2.rapidity()));
        histos["dijet_mass"]->fill((j1+j2).mass());
        histos["H_dijet_mass"]->fill((hmom+j1+j2).mass());
        histos["deltay_jj"]->fill(fabs(j1.rapidity()-j2.rapidity()));
        histos["deltay_H_jj"]->fill(fabs((hmom.rapidity()-(j1+j2).rapidity())));

        histos["loose"]->fill(0.);
        if (fabs(deltaPhi(hmom,j1+j2))>_dphiHjjtight) histos["tight"]->fill(0.);

        // introduce boolean whether to fill WBF and WBF2 observables
        bool wbf(fabs(j1.rapidity()-j2.rapidity())>_wbfdyjj &&
                 (j1+j2).mass()>_wbfmjj);

        bool wbf2(false);
        vector<Jet>::const_iterator itr1,itr2;
        for (itr1=RapJets.begin(); itr1!=RapJets.end()-1; ++itr1) {
          for (itr2=itr1+1; itr2!=RapJets.end(); ++itr2) {
            double dy_pair = fabs(itr1->momentum().rapidity()-itr2->momentum().rapidity());
            double m_pair = (itr1->momentum()+itr2->momentum()).mass();
            if (dy_pair>_wbfdyjj && m_pair>_wbfmjj) wbf2=true;
          }
        }

        // fill WBF histograms
        if (wbf) {
          histos["deltaphi2_VBF"]->fill(phi2);
          histos["H_jj_pT_VBF"]->fill(hmom.pT()/GeV);
          histos["deltaphi_jj_VBF"]->fill(deltaPhi(j1,j2));
          histos["deltaphi_Hjj_VBF"]->fill(deltaPhi(hmom,j1+j2));

          histos["NJet_excl_30_VBF"]->fill(PTJets.size());
          for (size_t i(0);i<4;++i) {
            if (PTJets.size()>=i) histos["NJet_incl_30_VBF"]->fill(i);
          }

          histos["loose"]->fill(1.);
          if (fabs(deltaPhi(hmom,j1+j2))>_dphiHjjtight) histos["tight"]->fill(1.);
        }

        // ** 15/5/2015 **

        // WBF2: look for any pair that fulfils the requirements

        // fill WBF2 histograms
        if (wbf2) {
          histos["deltaphi2_VBF2"]->fill(phi2);
          histos["H_jj_pT_VBF2"]->fill(hmom.pT()/GeV);
          histos["deltaphi_jj_VBF2"]->fill(deltaPhi(j1,j2));
          histos["deltaphi_Hjj_VBF2"]->fill(deltaPhi(hmom,j1+j2));

          histos["NJet_excl_30_VBF2"]->fill(PTJets.size());
          for (size_t i(0);i<4;++i) {
            if (PTJets.size()>=i) histos["NJet_incl_30_VBF2"]->fill(i);
          }

          histos["loose"]->fill(2.);
          if (fabs(deltaPhi(hmom,j1+j2))>_dphiHjjtight) histos["tight"]->fill(2.);
        }

        // END ** 15/5/2015 **

        double tauJet2 = sqrt(sqr(j2.pT()) + sqr(j2.mass()))/
          (2.0*cosh(j2.rapidity() - hmom.rapidity()));
        histos["tau_jet2"]->fill(tauJet2/GeV);

        // njets == 2;
        if (PTJets.size()==2) {
          histos["deltaphi_jj_excl"]->fill(deltaPhi(j1,j2));
          histos["deltaphi_Hjj_excl"]->fill(deltaPhi(hmom,j1+j2));
          histos["Hjj_pT_excl"]->fill((hmom+j1+j2).pT()/GeV);
          histos["H_jj_pT_excl"]->fill(hmom.pT()/GeV);
          histos["jet2_pT_excl"]->fill(j2.pT()/GeV);
          // again, why?
          // 6/2 added fill for j3_pT 0-30 GeV, i.e. no jet3
          // histos["jet3_pT_incl"]->fill(10);
        }


        // IVAN *******************************************************

        const size_t nj(PTJets.size());
        double jjpT_dy(0.), jjdy_dy(0.);

        jjpT_dy = fabs(j1.rapidity() - j2.rapidity());
        jjdy_dy = fabs(jf.rapidity() - jb.rapidity());

        // Fill histograms ----------------------

        // ** Jet pT for different dijet rapidity separations

        histos["jjpT_dy"]->fill(jjpT_dy);
        if      (nj==2) histos["jjpT_dy_2j_excl"]->fill(jjpT_dy);
        else if (nj==3) histos["jjpT_dy_3j_excl"]->fill(jjpT_dy);

        histos["jjdy_dy"]->fill(jjdy_dy);
        if      (nj==2) histos["jjdy_dy_2j_excl"]->fill(jjdy_dy);
        else if (nj==3) histos["jjdy_dy_3j_excl"]->fill(jjdy_dy);

        for (size_t j(0), nj(min(PTJets.size(),3lu)); j<nj; ++j) {
          const double pT = PTJets[j].momentum().pT();

          for (int y=1; y<=6; ++y) {
            stringstream ss;
            ss << "jet" << j+1 << "_pT_jjpT_mindy" << y;
            if ( jjpT_dy > (j+1) ) histos[ss.str()]->fill(pT);
            else break;
          }

          for (int y=1; y<=6; ++y) {
            stringstream ss;
            ss << "jet" << j+1 << "_pT_jjdy_mindy" << y;
            if ( jjdy_dy > (j+1) ) histos[ss.str()]->fill(pT);
            else break;
          }
        }


        // ** jet veto

        double ydists(100.);
        double ycenter=(jb.rapidity()+jf.rapidity())/2.;
        double ydistt;

        vector<Jet>::const_iterator jitr;
        for (jitr=RapJets.begin()+1; jitr!=RapJets.end()-1; ++jitr) {
          ydistt=fabs(jitr->momentum().rapidity()-ycenter);
          if (ydistt<ydists) ydists=ydistt;
        }

        // ydists is now the smallest distance between the centre of the
        // tagging jets and any possible further jet
        // (100 in case of no further jets)
        // Now fill the jet veto histogram
        Histo1DPtr const hcentraljveto = histos["xs_central_jet_veto"];
        Histo1DPtr const hcentraljveto_VBF = histos["xs_central_jet_veto_VBF"];
        Histo1DPtr const hcentraljveto_VBF2 = histos["xs_central_jet_veto_VBF2"];

        for (int i=0, n=hcentraljveto->numBins(); i<n; ++i) {
          if (hcentraljveto->bin(i).xMin() < ydists) {
            double xwidth=hcentraljveto->bin(i).xWidth();
            hcentraljveto->fillBin(i, xwidth);

            if (wbf)  hcentraljveto_VBF->fillBin(i, xwidth);
            // ** 15/5/2015 **
            if (wbf2) hcentraljveto_VBF2->fillBin(i, xwidth);
            // END ** 15/5/2015 **
          }
        }

        // rapidity distance between forward and backward jets
        histos["jjfb_dy"]->fill(fabs(jf.rapidity()-jb.rapidity()));

        // END IVAN ***************************************************
      }

      // njets > 2;
      if (PTJets.size()>2) {
        const FourMomentum& j3(PTJets[2].momentum());

        histos["H_jjj_pT_incl"]->fill(hmom.pT()/GeV);
        histos["jet3_pT_incl"]->fill(j3.pT()/GeV);
        histos["jet3_y"]->fill(fabs(j3.rapidity()));

        double tauJet3 = sqrt(sqr(j3.pT()) + sqr(j3.mass()))/
          (2.0*cosh(j3.rapidity() - hmom.rapidity()));
        histos["tau_jet3"]->fill(tauJet3/GeV);

        if (PTJets.size()==3) {
          histos["H_jjj_pT_excl"]->fill(hmom.pT()/GeV);
          histos["jet3_pT_excl"]->fill(j3.pT()/GeV);
        }
      }

      double HT_jets(0.),HT_all(0.);
      for(size_t i(0);i<PTJets.size();i++) HT_jets += PTJets[i].momentum().pT();
      HT_all=HT_jets+hmom.Et();
      histos["HT_jets"]->fill(HT_jets);
      histos["HT_all"]->fill(HT_all);

      double max_tj=0;
      double sum_tj=0;
      for (const Jet& j : PTJets) {
        double tauJet = sqrt(sqr(j.momentum().pT()) + sqr(j.momentum().mass()))/
          (2.0*cosh(j.momentum().rapidity() - hmom.rapidity()));
        if (tauJet > _taujcut) {
          sum_tj+=tauJet;
          if (tauJet > max_tj) max_tj=tauJet;
        }
      }
      if (max_tj>0.) histos["tau_jet_max"]->fill(max_tj);
      if (sum_tj>0.) histos["sum_tau_jet"]->fill(sum_tj);

      if (PTJets.size()>=0) {
        if (PTJets.size()>=1) {
          const FourMomentum& j1(PTJets[0].momentum());
          if (j1.pT()>50.*GeV) {
            histos["NJet_incl_50"]->fill(1);
          }
          else {
            histos["NJet_excl_50"]->fill(0);
          }
        }
        if (PTJets.size()==1) {
          const FourMomentum& j1(PTJets[0].momentum());
          if (j1.pT()>50.*GeV) {
            histos["NJet_excl_50"]->fill(1);
          }
        }
        if (PTJets.size()>=2) {
          const FourMomentum& j2(PTJets[1].momentum());
          if (j2.pT()>50.*GeV) {
            histos["NJet_incl_50"]->fill(2);
          }
        }
        if (PTJets.size()==2) {
          const FourMomentum& j2(PTJets[1].momentum());
          if(j2.pT()>50.*GeV) {
            histos["NJet_excl_50"]->fill(2);
          }
        }
        if (PTJets.size()>=3) {
          const FourMomentum& j3(PTJets[2].momentum());
          if (j3.pT()>50.*GeV) {
            histos["NJet_incl_50"]->fill(3);
          }
        }
        if (PTJets.size()==3) {
          const FourMomentum& j3(PTJets[2].momentum());
          if (j3.pT()>50.*GeV) {
            histos["NJet_excl_50"]->fill(3);
          }
        }
      }
      if (PTJets.size()==0) {
        histos["NJet_incl_50"]->fill(0);
        histos["NJet_excl_50"]->fill(0);
      }
    }


    /// Finalize
    void finalize() {
      const double scalefactor = crossSection()/sumOfWeights();
      for (auto& key_hist : histos) scale(key_hist.second, scalefactor);
    }

  };


  DECLARE_RIVET_PLUGIN(MC_HJETS_LH15);

}
