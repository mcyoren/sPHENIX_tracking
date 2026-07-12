// FitRDCA.C
//
// Run with defaults:
//   root -l -b -q 'FitRDCA.C()'
//
// Or specify input file:
//   root -l -b -q 'FitRDCA.C("output/testpp_pcaz10.root")'
//
// Or specify input file and output prefix:
//   root -l -b -q \
//   'FitRDCA.C("output/testpp_pcaz10.root","output/testpp_pcaz10_rdca_fit")'
//
// Default outputs:
//   output/rdca_sector_fit.root
//   output/rdca_sector_fit.pdf

#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TF1.h>
#include <TFitResultPtr.h>
#include <TH2.h>
#include <TMath.h>
#include <TProfile.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <vector>

namespace
{
  constexpr int NSectors = 12;
  constexpr int NFinalParameters = 3 + NSectors;

  // Final fit parameters:
  //
  // p[0]       = main sine amplitude
  // p[1]       = constant offset
  // p[2]       = main sine phase
  // p[3]...p[14] = amplitudes of sin(12*phi) in sectors 0...11
  //
  // Sector 0 is centered at phi = 0.
  Double_t IndividualSectorFit(Double_t* xx, Double_t* p)
  {
    const double phi = xx[0];

    const double sectorWidth =
        2.0 * TMath::Pi() / static_cast<double>(NSectors);

    int sector = static_cast<int>(
        std::floor((phi + 0.5 * sectorWidth) / sectorWidth));

    sector %= NSectors;

    if (sector < 0)
    {
      sector += NSectors;
    }

    const double mainModulation =
        p[0] * std::sin(phi + p[2]);

    const double sectorModulation =
        p[3 + sector] * std::sin(12.0 * phi);

    return mainModulation + p[1] + sectorModulation;
  }

  // Set an allowed interval around the previous fitted value.
  //
  // For example, with relativeDeviation = 0.10:
  //
  // previous value = 0.20
  // new limits     = 0.18 to 0.22
  //
  // The minimum half-width prevents parameters close to zero from
  // becoming effectively fixed.
  void SetRelativeLimits(
      TF1* function,
      const int parameter,
      const double centerValue,
      const double relativeDeviation = 0.10,
      const double minimumHalfWidth = 1.0e-4)
  {
    double halfWidth =
        relativeDeviation * std::abs(centerValue);

    halfWidth =
        std::max(halfWidth, minimumHalfWidth);

    const double lower =
        centerValue - halfWidth;

    const double upper =
        centerValue + halfWidth;

    function->SetParameter(
        parameter,
        centerValue);

    function->SetParLimits(
        parameter,
        lower,
        upper);
  }

  // Obtain the TProfile.
  //
  // First try an already existing profile:
  //
  //   h_rDCAzero_vs_atan2_px_py_pfx
  //   h_rDCA_vs_atan2_px_py_pfx
  //
  // Otherwise create ProfileX from:
  //
  //   h_rDCAzero_vs_atan2_px_py
  //   h_rDCA_vs_atan2_px_py
  TProfile* GetProfile(
      TDirectory* directory,
      const TString& outputName)
  {
    if (!directory)
    {
      return nullptr;
    }

    const std::vector<TString> profileNames = {
        "h_rDCAzero_vs_atan2_px_py_pfx",
        "h_rDCA_vs_atan2_px_py_pfx"};

    for (const TString& name : profileNames)
    {
      TProfile* sourceProfile = nullptr;

      directory->GetObject(
          name,
          sourceProfile);

      if (!sourceProfile)
      {
        continue;
      }

      TProfile* clonedProfile =
          dynamic_cast<TProfile*>(
              sourceProfile->Clone(outputName));

      if (clonedProfile)
      {
        clonedProfile->SetDirectory(nullptr);
        return clonedProfile;
      }
    }

    const std::vector<TString> histogramNames = {
        "h_rDCAzero_vs_atan2_px_py",
        "h_rDCA_vs_atan2_px_py"};

    for (const TString& name : histogramNames)
    {
      TH2* histogram = nullptr;

      directory->GetObject(
          name,
          histogram);

      if (!histogram)
      {
        continue;
      }

      TProfile* profile =
          histogram->ProfileX(outputName);

      if (profile)
      {
        profile->SetDirectory(nullptr);
        return profile;
      }
    }

    return nullptr;
  }

  bool FitCategory(
      TFile* inputFile,
      TFile* outputFile,
      const TString& category,
      const TString& pdfName,
      const bool firstPage,
      const bool lastPage)
  {
    TDirectory* directory =
        inputFile->GetDirectory(category);

    if (!directory)
    {
      std::cerr
          << "Directory not found: "
          << category
          << std::endl;

      return false;
    }

    std::unique_ptr<TProfile> profile(
        GetProfile(
            directory,
            Form("profile_%s", category.Data())));

    if (!profile)
    {
      std::cerr
          << "Could not find rDCA profile or TH2 in directory: "
          << category
          << std::endl;

      return false;
    }

    const double xmin = -TMath::Pi();
    const double xmax =  TMath::Pi();

    profile->SetTitle(
        Form(
            "rDCA_{zero} vs atan2(p_{x},p_{y}) [%s];"
            "atan2(p_{x},p_{y});"
            "#LT rDCA_{zero} #GT",
            category.Data()));

    profile->SetMarkerStyle(20);
    profile->SetMarkerSize(0.65);
    profile->SetLineWidth(1);

    std::cout
        << "\n====================================================\n"
        << "Fitting category: " << category << "\n"
        << "===================================================="
        << std::endl;

    // ============================================================
    // Stage 1:
    //
    // [0]*sin(x)
    // ============================================================

    TF1 fit1(
        Form("fit1_%s", category.Data()),
        "[0]*sin(x)",
        xmin,
        xmax);

    fit1.SetParName(
        0,
        "Main amplitude");

    fit1.SetParameter(
        0,
        0.3);

    fit1.SetNpx(50000);

    std::cout
        << "Stage 1: [0]*sin(x)"
        << std::endl;

    profile->Fit(
        &fit1,
        "RQ0");

    // ============================================================
    // Stage 2:
    //
    // [0]*sin(x) + [1]
    // ============================================================

    TF1 fit2(
        Form("fit2_%s", category.Data()),
        "[0]*sin(x) + [1]",
        xmin,
        xmax);

    fit2.SetParName(
        0,
        "Main amplitude");

    fit2.SetParName(
        1,
        "Offset");

    fit2.SetParameter(
        0,
        fit1.GetParameter(0));

    fit2.SetParameter(
        1,
        profile->GetMean(2));

    fit2.SetNpx(50000);

    std::cout
        << "Stage 2: [0]*sin(x) + [1]"
        << std::endl;

    profile->Fit(
        &fit2,
        "RQ0");

    // ============================================================
    // Stage 3:
    //
    // [0]*sin(x+[2]) + [1]
    // ============================================================

    TF1 fit3(
        Form("fit3_%s", category.Data()),
        "[0]*sin(x+[2]) + [1]",
        xmin,
        xmax);

    fit3.SetParName(
        0,
        "Main amplitude");

    fit3.SetParName(
        1,
        "Offset");

    fit3.SetParName(
        2,
        "Main phase");

    fit3.SetParameter(
        0,
        fit2.GetParameter(0));

    fit3.SetParameter(
        1,
        fit2.GetParameter(1));

    fit3.SetParameter(
        2,
        0.0);

    fit3.SetParLimits(
        2,
        -TMath::Pi(),
        TMath::Pi());

    fit3.SetNpx(50000);

    std::cout
        << "Stage 3: [0]*sin(x+[2]) + [1]"
        << std::endl;

    for (int iteration = 0;
         iteration < 5;
         ++iteration)
    {
      std::cout
          << "  Stage-3 iteration "
          << iteration + 1
          << "/5"
          << std::endl;

      profile->Fit(
          &fit3,
          "RQ0");
    }

    // ============================================================
    // Stage 4:
    //
    // [0]*sin(x+[2]) + [1] + [3]*sin(12*x)
    // ============================================================

    TF1 commonFit(
        Form("commonFit_%s", category.Data()),
        "[0]*sin(x+[2]) + [1] + [3]*sin(12*x)",
        xmin,
        xmax);

    commonFit.SetParName(
        0,
        "Main amplitude");

    commonFit.SetParName(
        1,
        "Offset");

    commonFit.SetParName(
        2,
        "Main phase");

    commonFit.SetParName(
        3,
        "Common sector amplitude");

    commonFit.SetParameter(
        0,
        fit3.GetParameter(0));

    commonFit.SetParameter(
        1,
        fit3.GetParameter(1));

    commonFit.SetParameter(
        2,
        fit3.GetParameter(2));

    commonFit.SetParameter(
        3,
        0.1);

    commonFit.SetParLimits(
        2,
        -TMath::Pi(),
        TMath::Pi());

    commonFit.SetNpx(50000);

    std::cout
        << "Stage 4: add common [3]*sin(12*x)"
        << std::endl;

    for (int iteration = 0;
         iteration < 10;
         ++iteration)
    {
      std::cout
          << "  Common-modulation iteration "
          << iteration + 1
          << "/10"
          << std::endl;

      profile->Fit(
          &commonFit,
          "RQ0");
    }

    const double initialMainAmplitude =
        commonFit.GetParameter(0);

    const double initialOffset =
        commonFit.GetParameter(1);

    const double initialMainPhase =
        commonFit.GetParameter(2);

    const double initialSectorAmplitude =
        commonFit.GetParameter(3);

    std::cout
        << "\nStabilized common-fit parameters:\n"
        << "  Main amplitude          = "
        << initialMainAmplitude
        << "\n  Offset                  = "
        << initialOffset
        << "\n  Main phase              = "
        << initialMainPhase
        << "\n  Common sector amplitude = "
        << initialSectorAmplitude
        << std::endl;

    // ============================================================
    // Stage 5:
    //
    // Replace the common sin(12*x) amplitude with 12 separate
    // sector amplitudes.
    //
    // In the first iteration, every sector amplitude may vary by
    // +/-10% around the common fitted amplitude.
    //
    // Before every later iteration, the allowed interval is moved
    // to +/-10% around that sector's result from the previous fit.
    // ============================================================

    TF1* finalFit = new TF1(
        Form("finalFit_%s", category.Data()),
        IndividualSectorFit,
        xmin,
        xmax,
        NFinalParameters);

    finalFit->SetParName(
        0,
        "Main amplitude");

    finalFit->SetParName(
        1,
        "Offset");

    finalFit->SetParName(
        2,
        "Main phase");

    finalFit->SetParameter(
        0,
        initialMainAmplitude);

    finalFit->SetParameter(
        1,
        initialOffset);

    finalFit->SetParameter(
        2,
        initialMainPhase);

    finalFit->SetParLimits(
        2,
        -TMath::Pi(),
        TMath::Pi());

    for (int sector = 0;
         sector < NSectors;
         ++sector)
    {
      const int parameter =
          3 + sector;

      finalFit->SetParName(
          parameter,
          Form("Sector %02d amplitude", sector));

      SetRelativeLimits(
          finalFit,
          parameter,
          initialSectorAmplitude,
          0.10);
    }

    // Large Npx is useful because the final function is piecewise.
    finalFit->SetNpx(100000);

    constexpr int NFinalIterations = 10;

    std::cout
        << "\nStage 5: 12 sector-dependent amplitudes\n"
        << "Each iteration allows each amplitude to move by +/-10%\n"
        << "relative to its result from the preceding iteration."
        << std::endl;

    for (int iteration = 0;
         iteration < NFinalIterations;
         ++iteration)
    {
      std::cout
          << "\n  Final-stage iteration "
          << iteration + 1
          << "/"
          << NFinalIterations
          << std::endl;

      if (iteration > 0)
      {
        for (int sector = 0;
             sector < NSectors;
             ++sector)
        {
          const int parameter =
              3 + sector;

          const double previousValue =
              finalFit->GetParameter(parameter);

          SetRelativeLimits(
              finalFit,
              parameter,
              previousValue,
              0.10);
        }
      }

      const TString fitOptions =
          (iteration == NFinalIterations - 1)
              ? "RSME"
              : "RQ0S";

      TFitResultPtr result =
          profile->Fit(
              finalFit,
              fitOptions);

      std::cout
          << "    Fit status = "
          << static_cast<int>(result)
          << "\n    chi2/NDF   = "
          << finalFit->GetChisquare()
          << " / "
          << finalFit->GetNDF()
          << std::endl;

      for (int sector = 0;
           sector < NSectors;
           ++sector)
      {
        const int parameter =
            3 + sector;

        std::cout
            << "    Sector "
            << sector
            << " amplitude = "
            << finalFit->GetParameter(parameter)
            << " +/- "
            << finalFit->GetParError(parameter)
            << std::endl;
      }
    }

    finalFit->SetLineColor(kRed);
    finalFit->SetLineWidth(2);
    finalFit->SetNpx(100000);

    std::cout
        << "\nFinal result for "
        << category
        << "\n"
        << "Chi2/NDF       = "
        << finalFit->GetChisquare()
        << " / "
        << finalFit->GetNDF()
        << "\n"
        << "Main amplitude = "
        << finalFit->GetParameter(0)
        << " +/- "
        << finalFit->GetParError(0)
        << "\n"
        << "Offset         = "
        << finalFit->GetParameter(1)
        << " +/- "
        << finalFit->GetParError(1)
        << "\n"
        << "Main phase     = "
        << finalFit->GetParameter(2)
        << " +/- "
        << finalFit->GetParError(2)
        << std::endl;

    for (int sector = 0;
         sector < NSectors;
         ++sector)
    {
      const int parameter =
          3 + sector;

      std::cout
          << "Sector "
          << sector
          << " amplitude = "
          << finalFit->GetParameter(parameter)
          << " +/- "
          << finalFit->GetParError(parameter)
          << std::endl;
    }

    // ============================================================
    // Draw result
    // ============================================================

    TCanvas canvas(
        Form("canvas_%s", category.Data()),
        category,
        1600,
        1000);

    canvas.SetLeftMargin(0.11);
    canvas.SetRightMargin(0.22);
    canvas.SetBottomMargin(0.12);
    canvas.SetTopMargin(0.10);

    profile->Draw("E1");
    finalFit->Draw("same");

    canvas.Modified();
    canvas.Update();

    TString pdfOutput =
        pdfName;

    if (firstPage)
    {
      pdfOutput += "(";
    }
    else if (lastPage)
    {
      pdfOutput += ")";
    }

    canvas.Print(pdfOutput);

    // ============================================================
    // Save profile, fit, and canvas
    // ============================================================

    outputFile->cd();

    TDirectory* outputDirectory =
        outputFile->GetDirectory(category);

    if (!outputDirectory)
    {
      outputDirectory =
          outputFile->mkdir(category);
    }

    outputDirectory->cd();

    profile->Write(
        profile->GetName(),
        TObject::kOverwrite);

    finalFit->Write(
        finalFit->GetName(),
        TObject::kOverwrite);

    canvas.Write(
        canvas.GetName(),
        TObject::kOverwrite);

    outputFile->cd();

    return true;
  }
}

void FitRDCA(
    const char* inputFileName = "output/testpp_pcaz10.root",
    const char* outputPrefix = "outputfit/rdca_sector_fit")
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);

  // Create output directory if it does not exist.
  gSystem->mkdir(
      "output",
      true);

  std::unique_ptr<TFile> inputFile(
      TFile::Open(
          inputFileName,
          "READ"));

  if (!inputFile || inputFile->IsZombie())
  {
    std::cerr
        << "Could not open input file: "
        << inputFileName
        << std::endl;

    return;
  }

  const TString outputRootName =
      Form("%s.root", outputPrefix);

  const TString outputPdfName =
      Form("%s.pdf", outputPrefix);

  std::unique_ptr<TFile> outputFile(
      TFile::Open(
          outputRootName,
          "RECREATE"));

  if (!outputFile || outputFile->IsZombie())
  {
    std::cerr
        << "Could not create output ROOT file: "
        << outputRootName
        << std::endl;

    return;
  }

  const std::vector<TString> categories = {
      "side0_qplus",
      "side0_qminus",
      "side1_qplus",
      "side1_qminus"};

  int successfulFits = 0;

  for (std::size_t index = 0;
       index < categories.size();
       ++index)
  {
    const bool success =
        FitCategory(
            inputFile.get(),
            outputFile.get(),
            categories[index],
            outputPdfName,
            index == 0,
            index == categories.size() - 1);

    if (success)
    {
      ++successfulFits;
    }
  }

  outputFile->Write();
  outputFile->Close();
  inputFile->Close();

  std::cout
      << "\n====================================================\n"
      << "Finished "
      << successfulFits
      << " of "
      << categories.size()
      << " fits.\n"
      << "Output ROOT file: "
      << outputRootName
      << "\n"
      << "Output PDF:       "
      << outputPdfName
      << "\n"
      << "===================================================="
      << std::endl;
}