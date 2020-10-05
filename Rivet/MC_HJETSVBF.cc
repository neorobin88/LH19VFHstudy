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
      for (size_t i = 0; i < DRS.size(); ++i) {
        FastJets fj(fs, FastJets::ANTIKT, DRS[i]);
        string dr = "_xxx";
        sprintf(&dr[1], "%0.1f", DRS[i]);
        declare(fj, "Jets"+dr);
      }
      PromptFinalState pfs(Cuts::abseta < 5);
      PromptFinalState higgses(Cuts::pid == PID::HIGGS && Cuts::absrap < 2.4);
      declare(pfs, "PromptFS");
      declare(higgses, "Higgses");


      // Histograms
      /// @todo Add log-scale histo copies?
      const doubles edges_njets = linspace(6, -0.5, 5.5); //  from 0 to 5, exclusive and inclusive
      const doubles edges_delta_y_jj = linspace(200, 0.0, 10.0); // 0.05 bins from 0-10 (can be recombined later if needed) (log and linear)
      const doubles edges_delta_phi_jj = linspace(30, 0.0, 1.0); // 0.1 bins from 0-pi/pi
      const doubles edges_m_jj = linspace(40, 0.0, 2000.0); // 50 GeV bins from 0-2000 GeV (log and linear)
      const doubles edges_ptn_ptm = linspace(20, 0.0, 1.0); // 0.05 bins from 0-1
      const doubles edges_ht = linspace(30, 0.0, 3000.0); // 100 GeV bins from 0-3 TeV (log and linear)
      const doubles edges_xs = linspace(10, 0.0, 0.5); // 0.05 bins from 0 to 0.5
      const doubles edges_pth = linspace(50, 0.0, 500.0); // 10 GeV bins from 0-500 GeV (log and linear)
      const doubles edges_pthj = linspace(40, 0.0, 200.0); // -pTHj; 5 GeV bins from 0-200 GeV (log and linear)
      const doubles edges_pthjj = linspace(20, 0.0, 100.0); // -pTHjj: 5 GeV bins from 0-100 GeV (log and linear)

      // Inclusive histograms (dR = 0.4 only)
      // [njets, delta_y_jj, m_jj,delta_phi_jj, HT, pTH, pTHj, pTHjj]
      string pre = "incl_";
      book(_h_incl["njets"], pre+"njets", edges_njets);
      book(_h_incl["delta_y_jj12"], pre+"delta_y_jj12", edges_delta_y_jj);
      book(_h_incl["delta_phi_jj12"], pre+"delta_phi_jj12", edges_delta_phi_jj);
      book(_h_incl["m_jj12"], pre+"m_jj12", edges_m_jj);
      book(_h_incl["ht"], pre+"HT", edges_ht);
      book(_h_incl["pth"], pre+"pTH", edges_pth);
      book(_h_incl["pthj1"], pre+"pTHj1", edges_pthj);
      book(_h_incl["pthjj12"], pre+"pTHjj12", edges_pthjj);

      for (size_t ir = 0; ir < DRS.size(); ++ir) {
        string dr = "xxx";
        sprintf(&dr[0], "%0.1f", DRS[ir]);

        // Per-dR, per-pTH, per-dy histograms
        // [njets, delta_y_jj12, m_jj12, delta_phi_jj12, delta_y_jjfb,
        //  pT2/pT1, pT3/pT1, xH, x1, x2, x3, pTH, pTHj, pTHjj]
        for (size_t ih = 0; ih < PTHCUTS.size(); ++ih) {
          const string pth = PTHCUTS[ih];
          for (size_t iy = 0; iy < DY12CUTS.size(); ++iy) {
            const string dy = DY12CUTS[iy];
            const string pre = "drptdy_" + dr + "_" + pth + "_" + dy + "_";
            const string suff = "";
            book(_h_dr_pth_dy[ir][ih][iy]["njets"], pre+"njets"+suff, edges_njets);
            book(_h_dr_pth_dy[ir][ih][iy]["delta_y_jj12"], pre+"delta_y_jj12"+suff, edges_delta_y_jj);
            book(_h_dr_pth_dy[ir][ih][iy]["delta_phi_jj12"], pre+"delta_phi_jj12"+suff, edges_delta_phi_jj);
            book(_h_dr_pth_dy[ir][ih][iy]["m_jj12"], pre+"m_jj12"+suff, edges_m_jj);
            book(_h_dr_pth_dy[ir][ih][iy]["delta_y_jjfb"], pre+"delta_y_jjfb"+suff, edges_delta_y_jj);
            book(_h_dr_pth_dy[ir][ih][iy]["pt2_pt1"], pre+"pt2_pt1"+suff, edges_ptn_ptm);
            book(_h_dr_pth_dy[ir][ih][iy]["pt3_pt1"], pre+"pt3_pt1"+suff, edges_ptn_ptm);
            book(_h_dr_pth_dy[ir][ih][iy]["xh"], pre+"xh"+suff, edges_xs);
            book(_h_dr_pth_dy[ir][ih][iy]["x1"], pre+"x1"+suff, edges_xs);
            book(_h_dr_pth_dy[ir][ih][iy]["x2"], pre+"x2"+suff, edges_xs);
            book(_h_dr_pth_dy[ir][ih][iy]["x3"], pre+"x3"+suff, edges_xs);
            book(_h_dr_pth_dy[ir][ih][iy]["ht"], pre+"ht"+suff, edges_xs);
            book(_h_dr_pth_dy[ir][ih][iy]["pth"], pre+"pth"+suff, edges_pth);
            book(_h_dr_pth_dy[ir][ih][iy]["pthj1"], pre+"pthj1"+suff, edges_pthj);
            book(_h_dr_pth_dy[ir][ih][iy]["pthjj12"], pre+"pthjj12"+suff, edges_pthjj);
          }
        }

        // Per-dR, per-rescuts histograms
        // [delta_y_jj12, m_jj12, delta_phi_jj12, delta_y_jjfb]
        for (size_t iv = 0; iv < 2; ++iv) {
          const string res = RESNAMES[iv];
          const string pre = "drres_" + dr + "_" + pth + "_" + res + "_";
          const string suff = "";
          book(_h_dr_res[ir][ih][iy]["delta_y_jj12"], pre+"delta_y_jj12"+suff, edges_delta_y_jj);
          book(_h_dr_res[ir][ih][iy]["delta_phi_jj12"], pre+"delta_phi_jj12"+suff, edges_delta_phi_jj);
          book(_h_dr_res[ir][ih][iy]["m_jj12"], pre+"m_jj12"+suff, edges_m_jj);
          book(_h_dr_res[ir][ih][iy]["delta_y_jjfb"], pre+"delta_y_jjfb"+suff, edges_delta_y_jj);
        }

        // Higgs pT histogram with ATLAS VBF cuts
        pre = "atlas_";
        book(_h_atlas["pth"], pre+"pth", edges_pth);

      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Get Higgs
      const Particles higgses = apply<ParticleFinder>(event, "Higgses").particles();
      if (higgses.empty()) continue;
      if (higgses.size() > 1) continue;
      const Particle higgs = higgses[0];
      const double ptH = higgs.pT();
      const double yh  = higgs.rap();

      for (size_t ir = 0; ir < DRS.size(); ++j) {
        string dr = "xxx";
        sprintf(&dr[0], "%0.1f", DRS[ir]);

        // Get jets, require dijet, and compute HTs
        const Jets jets = apply<FastJets>(event, "Jets"+dr)
          .jetsByPt(Cuts::pT > 30*GeV && Cuts::absrap < 4.4);
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
        const double dR12   = add_quad(dy12, dphi12);
        const double y1  = jets[0].rap();
        const double y2  = jets[1].rap();

        // Get most fwd/bwd jet pair
        const Jets jetsbyrap = sortBy(jets, rapGtr); //< sorts from fwd to bwd
        const Jet& jf = jets.front();
        const Jet& jb = jets.back();
        const FourMomentum dijetfb = jf.mom() + jb.mom();
        const double mfb = dijetfb.mass();
        const double ptfb = dijetfb.pT();
        const double yfb = dijetfb.rap();
        const double dyfb = std::abs(jf.rap() - jb.rap());
        const double dphifb = deltaPhi(jf, jb);
        const double dRfb   = add_quad(dyfb, dphifb);
        const double y1  = jf.rap();
        const double y2  = jb.rap();

        // Get leading jet + H
        const FourMomentum jH = jets[0].mom() + higgs.mom();
        const double ptjH = jH.pT();

        // Get leading dijet + H
        const FourMomentum jjH = dijet + higgs.mom();
        const double ptjjH = jjH.pT();
        const double dphijjH = deltaPhi(higgs, dijet);

        // Get VH resonance-cut status
        const double minmasswindow = 50*GeV, maxmasswindow = 150*GeV;
        bool novh = !inRange(m12, minmasswindow, maxmasswindow);
        if (novh) {
          const double m1 = jets[0].mom().mass();
          const double m2 = jets[1].mom().mass();
          if (m1 > minmasswindow && m1 < maxmasswindow) novh = false;
          if (m2 > minmasswindow && m2 < maxmasswindow) novh = false;

          if (njets > 2) {
            const double m3 = jets[2].mom().mass();
            const double m13 = (jets[0].mom() + jets[2].mom()).mass();
            const double m23 = (jets[1].mom() + jets[2].mom()).mass();
            const double m123 = (jets[0].mom() +jets[1].mom() + jets[2].mom()).mass();
            if (m3 > minmasswindow && m3 < maxmasswindow) novh = false;
            if (m13 > minmasswindow && m13 < maxmasswindow) novh = false;
            if (m23 > minmasswindow && m23 < maxmasswindow) novh = false;
            if (m123 > minmasswindow && m123 < maxmasswindow) novh = false;
          }
        }


        // Inclusive plots, for dR = 0.4 only
        /// @todo Aren't these duplicated in the cut-combination plots?
        if (ir == 0) { //< dR = 0.4
          _h_incl["njets"]->fill(njets);
          _h_incl["delta_y_jj12"]->fill(dy12);
          _h_incl["delta_phi_jj12"]->fill(dphi12/M_PI);
          _h_incl["m_jj12"]->fill(m12/GeV);
          _h_incl["ht"]->fill(ht/GeV);
          _h_incl["pth"]->fill(ptH/GeV);
          _h_incl["pthj1"]->fill(ptjH/GeV);
          _h_incl["pthjj12"]->fill(ptjjH/GeV);
        }

        // Standard histograms for pTH and dy12 cuts
        for (size_t ih = 0; ih < PTHS.size(); ++ih) {
          if (pTH < PTHCUTS[ih]*GeV) continue;
          for (size_t iy = 0; iy < 2; ++iy) {
            if (dy12 < DY12CUTS[iy]) continue;
            // _c_xs[j][i]->fill();
            _h_incl[ir][ih][iy]["njets"]->fill(njets);
            _h_incl[ir][ih][iy]["delta_y_jj12"]->fill(dy12);
            _h_incl[ir][ih][iy]["delta_phi_jj12"]->fill(dphi12/M_PI);
            _h_incl[ir][ih][iy]["m_jj12"]->fill(m12/GeV);
            _h_incl[ir][ih][iy]["delta_y_jjfb"]->fill(dyfb);
            _h_incl[ir][ih][iy]["pt2_pt1"]->fill(ptjjH/GeV);
            _h_incl[ir][ih][iy]["pt3_pt1"]->fill(ptjjH/GeV);
            _h_incl[ir][ih][iy]["ht"]->fill(ht/GeV);
            _h_incl[ir][ih][iy]["pth"]->fill(ptH/GeV);
            _h_incl[ir][ih][iy]["pthj1"]->fill(ptjH/GeV);
            _h_incl[ir][ih][iy]["pthjj12"]->fill(ptjjH/GeV);
            _h_incl[ir][ih][iy]["xh"]->fill(ptH/ht);
            _h_incl[ir][ih][iy]["x1"]->fill(jets[0].pT()/ht);
            _h_incl[ir][ih][iy]["x2"]->fill(jets[1].pT()/ht);
            if (njets > 2)
              _h_incl[ir][ih][iy]["x3"]->fill(jets[2].pT()/ht);
          }
        }

        // Histograms with/without the VH resonance cut
        for (size_t iv = 0; iv < 2; ++iv) {
          if (bool(iv) == novh) continue; // if passed, 0 = !res, 1 = res
          _h_dr_res[ir][iv]["delta_y_jj12"]->fill(dy12);
          _h_dr_res[ir][iv]["delta_phi_jj12"]->fill(dphi12/M_PI);
          _h_dr_res[ir][iv]["m_jj12"]->fill(m12/GeV);
          _h_dr_res[ir][iv]["delta_y_jjfb"]->fill(dyfb);
        }

        // Final cuts used for ATLAS VBF measurements in a Higgs pT histogram
        //   delta_y_jj>3.5 (where the two jets are the two highest pT jets),
        //   m_jj > 600 GeV (again using the two highest pT jets) and
        //   delta_phi_H_jj>2.7.
        if (dy12 > 3.5 && m12 > 600*GeV && dphijjH > 2.7) {
          _h_atlas["pth"]->fill(ptH/GeV);
        }

      } // DR loop
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      const double sf = crossSection()/femtobarn/sumOfWeights();
      scale(_h_incl, sf);
      for (size_t ir = 0; ir < DRS.size(); ++ir) {
        for (size_t ih = 0; ih < PTHCUTS.size(); ++ih) {
          for (size_t iy = 0; iy < DY12CUTS.size(); ++iy) {
            scale(_h_dr_pth_dy[ir][ih][iy], sf);
          }
        }
        for (size_t iv = 0; iv < 2; ++iv) {
          scale(_h_dr_res[ir][iv], sf);
        }
      }
      scale(_h_atlas, sf);
    }


    // Histograms
    map<string, Histo1DPtr> _h_incl, _h_dr_pth_dy[3][3][2], _h_dr_res[3][2], _h_atlas;

    // Cut values for standard histogram sets (other than m12 and dy12)
    static const vector<double> DRS;
    static const vector<double> PTHCUTS;
    static const vector<double> DY12CUTS;
    static const vector<double> RESNAMES;

  };


  // Static const initializers
  const vector<double> MC_HJETSVBF::DRS = {0.4, 0.7, 1.0};
  const vector<double> MC_HJETSVBF::PTHCUTS = {0., 200., 500.};
  const vector<double> MC_HJETSVBF::DY12CUTS = {0.0, 1.0};
  const vector<string> MC_HJETSVBF::RESNAMES = {"nores", "res"};


  DECLARE_RIVET_PLUGIN(MC_HJETSVBF);

}
