// FitRDCASlices.C
//
// Run:
//   root -l -b -q 'FitRDCASlices.C()'
//
// Custom input:
//   root -l -b -q \
//   'FitRDCASlices.C("output/testpp_pcaz10.root")'
//
// Custom input and output prefix:
//   root -l -b -q \
//   'FitRDCASlices.C("output/testpp_pcaz10.root",
//                     "output/testpp_pcaz10_slice_fit")'
//
// Default outputs:
//   output/rdca_slice_fit.root
//   output/rdca_slice_fit.pdf

#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TH2.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
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
  constexpr int NFinalParameters = 5 + NSectors;

  // ------------------------------------------------------------
  // Configuration
  // ------------------------------------------------------------

  constexpr int PhiBinGroup = 2;

  constexpr int MinimumSliceEntries = 50;

  constexpr int NGaussianIterations = 4;

  constexpr int NCommonFitIterations = 10;

  constexpr int NFinalFitIterations = 10;

  constexpr double SliceFitNSigma = 2.0;

  constexpr double RelativeStep = 0.10;

  constexpr double MinimumParameterHalfWidth = 1.0e-4;

  constexpr double MinimumAcceptedSigma = 0.01;

  constexpr double MaximumAcceptedSigma = 2.0;

  // Set true to use a shared-mean double Gaussian after the
  // iterative single-Gaussian fit.
  constexpr bool UseDoubleGaussian = true;

  // ------------------------------------------------------------
  // Final piecewise function
  //
  // p[0]       = main sine amplitude
  // p[1]       = offset
  // p[2]       = main sine phase
  // p[3]       = common 12-fold phase
  // p[4]       = optional second-harmonic amplitude
  // p[5]...p[16] = sector-dependent 12-fold amplitudes
  //
  // Final model:
  //
  // A1 sin(phi + phase1)
  // + offset
  // + Asector sin(12 phi + phase12)
  // + A2 sin(2 phi)
  //
  // The second harmonic is fixed at zero by default below.
  // ------------------------------------------------------------

  Double_t FinalSectorFunction(Double_t* xx, Double_t* p)
  {
    const double phi = xx[0];

    const double sectorWidth =
        2.0 * TMath::Pi() / static_cast<double>(NSectors);

    // Sector zero is centered at phi = 0.
    int sector = static_cast<int>(
        std::floor((phi + 0.5 * sectorWidth) / sectorWidth));

    sector %= NSectors;

    if (sector < 0)
    {
      sector += NSectors;
    }

    const double main =
        p[0] * std::sin(phi + p[2]);

    const double sectorTerm =
        p[5 + sector] * std::sin(12.0 * phi + p[3]);

    const double secondHarmonic =
        p[4] * std::sin(2.0 * phi);

    return p[1] + main + sectorTerm + secondHarmonic;
  }

  // ------------------------------------------------------------
  // Shared-mean double Gaussian
  //
  // p[0] = narrow normalization
  // p[1] = shared mean
  // p[2] = narrow sigma
  // p[3] = wide normalization
  // p[4] = wide sigma
  // ------------------------------------------------------------

  Double_t SharedMeanDoubleGaussian(Double_t* xx, Double_t* p)
  {
    const double x = xx[0];

    const double narrowArgument =
        (x - p[1]) / p[2];

    const double wideArgument =
        (x - p[1]) / p[4];

    return
        p[0] * std::exp(-0.5 * narrowArgument * narrowArgument)
        +
        p[3] * std::exp(-0.5 * wideArgument * wideArgument);
  }

  // ------------------------------------------------------------
  // Set limits around the value from the previous iteration.
  // ------------------------------------------------------------

  void SetRelativeLimits(
      TF1* function,
      const int parameter,
      const double center,
      const double relativeWidth = RelativeStep)
  {
    double halfWidth =
        relativeWidth * std::abs(center);

    halfWidth =
        std::max(halfWidth, MinimumParameterHalfWidth);

    function->SetParameter(parameter, center);

    function->SetParLimits(
        parameter,
        center - halfWidth,
        center + halfWidth);
  }

  // ------------------------------------------------------------
  // Locate the 2D input histogram.
  // ------------------------------------------------------------

  TH2* GetInputHistogram(TDirectory* directory)
  {
    if (!directory)
    {
      return nullptr;
    }

    const std::vector<TString> names = {
        "h_rDCAzero_vs_atan2_px_py",
        "h_rDCA_vs_atan2_px_py"};

    for (const TString& name : names)
    {
      TH2* histogram = nullptr;

      directory->GetObject(name, histogram);

      if (histogram)
      {
        return histogram;
      }
    }

    return nullptr;
  }

  // ------------------------------------------------------------
  // Perform iterative single-Gaussian slice fit.
  // ------------------------------------------------------------

  bool FitSliceGaussian(
      TH1D* slice,
      double& fittedMean,
      double& fittedMeanError,
      double& fittedSigma)
  {
    if (!slice)
    {
      return false;
    }

    if (slice->GetEntries() < MinimumSliceEntries)
    {
      return false;
    }

    const int maximumBin =
        slice->GetMaximumBin();

    const double peakPosition =
        slice->GetXaxis()->GetBinCenter(maximumBin);

    double initialSigma =
        slice->GetRMS();

    if (!std::isfinite(initialSigma) ||
        initialSigma < MinimumAcceptedSigma)
    {
      initialSigma = 0.3;
    }

    initialSigma =
        std::min(initialSigma, 1.0);

    double fitMinimum =
        peakPosition - 1.5 * initialSigma;

    double fitMaximum =
        peakPosition + 1.5 * initialSigma;

    fitMinimum =
        std::max(
            fitMinimum,
            slice->GetXaxis()->GetXmin());

    fitMaximum =
        std::min(
            fitMaximum,
            slice->GetXaxis()->GetXmax());

    TF1 gaussian(
        Form("slice_gaussian_%s", slice->GetName()),
        "gaus",
        fitMinimum,
        fitMaximum);

    gaussian.SetParameters(
        slice->GetMaximum(),
        peakPosition,
        initialSigma);

    gaussian.SetParLimits(
        2,
        MinimumAcceptedSigma,
        MaximumAcceptedSigma);

    TFitResultPtr result =
        slice->Fit(&gaussian, "RQ0S");

    if (static_cast<int>(result) != 0)
    {
      return false;
    }

    for (int iteration = 0;
         iteration < NGaussianIterations;
         ++iteration)
    {
      const double mean =
          gaussian.GetParameter(1);

      const double sigma =
          std::abs(gaussian.GetParameter(2));

      if (!std::isfinite(mean) ||
          !std::isfinite(sigma) ||
          sigma < MinimumAcceptedSigma ||
          sigma > MaximumAcceptedSigma)
      {
        return false;
      }

      fitMinimum =
          mean - SliceFitNSigma * sigma;

      fitMaximum =
          mean + SliceFitNSigma * sigma;

      fitMinimum =
          std::max(
              fitMinimum,
              slice->GetXaxis()->GetXmin());

      fitMaximum =
          std::min(
              fitMaximum,
              slice->GetXaxis()->GetXmax());

      gaussian.SetRange(
          fitMinimum,
          fitMaximum);

      result =
          slice->Fit(&gaussian, "RQ0S");

      if (static_cast<int>(result) != 0)
      {
        return false;
      }
    }

    fittedMean =
        gaussian.GetParameter(1);

    fittedMeanError =
        gaussian.GetParError(1);

    fittedSigma =
        std::abs(gaussian.GetParameter(2));

    if (!std::isfinite(fittedMean) ||
        !std::isfinite(fittedMeanError) ||
        !std::isfinite(fittedSigma))
    {
      return false;
    }

    if (fittedMeanError <= 0.0 ||
        fittedSigma < MinimumAcceptedSigma ||
        fittedSigma > MaximumAcceptedSigma)
    {
      return false;
    }

    return true;
  }

  // ------------------------------------------------------------
  // Optional shared-mean double-Gaussian refinement.
  // ------------------------------------------------------------

  bool FitSliceDoubleGaussian(
      TH1D* slice,
      const double gaussianMean,
      const double gaussianMeanError,
      const double gaussianSigma,
      double& fittedMean,
      double& fittedMeanError)
  {
    if (!slice)
    {
      return false;
    }

    const double fitMinimum =
        std::max(
            gaussianMean - 4.0 * gaussianSigma,
            slice->GetXaxis()->GetXmin());

    const double fitMaximum =
        std::min(
            gaussianMean + 4.0 * gaussianSigma,
            slice->GetXaxis()->GetXmax());

    TF1 doubleGaussian(
        Form("double_gaussian_%s", slice->GetName()),
        SharedMeanDoubleGaussian,
        fitMinimum,
        fitMaximum,
        5);

    doubleGaussian.SetParNames(
        "Core norm",
        "Common mean",
        "Core sigma",
        "Tail norm",
        "Tail sigma");

    doubleGaussian.SetParameters(
        slice->GetMaximum(),
        gaussianMean,
        std::max(0.7 * gaussianSigma, MinimumAcceptedSigma),
        0.20 * slice->GetMaximum(),
        std::min(2.0 * gaussianSigma, MaximumAcceptedSigma));

    doubleGaussian.SetParLimits(
        0,
        0.0,
        10.0 * slice->GetMaximum());

    doubleGaussian.SetParLimits(
        1,
        gaussianMean - gaussianSigma,
        gaussianMean + gaussianSigma);

    doubleGaussian.SetParLimits(
        2,
        MinimumAcceptedSigma,
        std::max(gaussianSigma, 0.05));

    doubleGaussian.SetParLimits(
        3,
        0.0,
        10.0 * slice->GetMaximum());

    doubleGaussian.SetParLimits(
        4,
        std::max(gaussianSigma, 0.02),
        MaximumAcceptedSigma);

    TFitResultPtr result =
        slice->Fit(&doubleGaussian, "RQ0S");

    if (static_cast<int>(result) != 0)
    {
      return false;
    }

    // Repeat a few times using the preceding result.
    for (int iteration = 0;
         iteration < 3;
         ++iteration)
    {
      result =
          slice->Fit(&doubleGaussian, "RQ0S");

      if (static_cast<int>(result) != 0)
      {
        return false;
      }
    }

    fittedMean =
        doubleGaussian.GetParameter(1);

    fittedMeanError =
        doubleGaussian.GetParError(1);

    if (!std::isfinite(fittedMean) ||
        !std::isfinite(fittedMeanError) ||
        fittedMeanError <= 0.0)
    {
      fittedMean = gaussianMean;
      fittedMeanError = gaussianMeanError;
      return false;
    }

    return true;
  }

  // ------------------------------------------------------------
  // Build TGraphErrors from fitted Y slices.
  // ------------------------------------------------------------

  TGraphErrors* BuildSliceGraph(
      TH2* histogram,
      const TString& graphName,
      TDirectory* diagnosticDirectory)
  {
    if (!histogram)
    {
      return nullptr;
    }

    auto* graph =
        new TGraphErrors();

    graph->SetName(graphName);

    graph->SetTitle(
        Form(
            "%s;atan2(p_{x},p_{y});fitted rDCA center",
            histogram->GetTitle()));

    int graphPoint = 0;

    const int numberOfXBins =
        histogram->GetNbinsX();

    for (int firstBin = 1;
         firstBin <= numberOfXBins;
         firstBin += PhiBinGroup)
    {
      const int lastBin =
          std::min(
              firstBin + PhiBinGroup - 1,
              numberOfXBins);

      const TString sliceName =
          Form(
              "%s_slice_%03d_%03d",
              graphName.Data(),
              firstBin,
              lastBin);

      std::unique_ptr<TH1D> slice(
          histogram->ProjectionY(
              sliceName,
              firstBin,
              lastBin,
              "e"));

      if (!slice)
      {
        continue;
      }

      slice->SetDirectory(nullptr);

      if (slice->GetEntries() < MinimumSliceEntries)
      {
        continue;
      }

      const double phiLow =
          histogram->GetXaxis()->GetBinLowEdge(firstBin);

      const double phiHigh =
          histogram->GetXaxis()->GetBinUpEdge(lastBin);

      const double phi =
          0.5 * (phiLow + phiHigh);

      const double phiError =
          0.5 * (phiHigh - phiLow);

      double gaussianMean = 0.0;
      double gaussianMeanError = 0.0;
      double gaussianSigma = 0.0;

      const bool gaussianSuccess =
          FitSliceGaussian(
              slice.get(),
              gaussianMean,
              gaussianMeanError,
              gaussianSigma);

      if (!gaussianSuccess)
      {
        continue;
      }

      double finalMean =
          gaussianMean;

      double finalMeanError =
          gaussianMeanError;

      if (UseDoubleGaussian)
      {
        double doubleGaussianMean = gaussianMean;
        double doubleGaussianMeanError = gaussianMeanError;

        const bool doubleGaussianSuccess =
            FitSliceDoubleGaussian(
                slice.get(),
                gaussianMean,
                gaussianMeanError,
                gaussianSigma,
                doubleGaussianMean,
                doubleGaussianMeanError);

        if (doubleGaussianSuccess)
        {
          finalMean =
              doubleGaussianMean;

          finalMeanError =
              doubleGaussianMeanError;
        }
      }

      // Prevent unrealistically tiny errors from dominating the
      // global modulation fit.
      const double minimumGlobalError = 0.005;

      finalMeanError =
          std::max(
              finalMeanError,
              minimumGlobalError);

      graph->SetPoint(
          graphPoint,
          phi,
          finalMean);

      graph->SetPointError(
          graphPoint,
          phiError,
          finalMeanError);

      ++graphPoint;

      if (diagnosticDirectory)
      {
        diagnosticDirectory->cd();
        slice->Write();
      }
    }

    if (graph->GetN() == 0)
    {
      delete graph;
      return nullptr;
    }

    graph->SetMarkerStyle(20);
    graph->SetMarkerSize(0.8);
    graph->SetLineWidth(1);

    return graph;
  }

  // ------------------------------------------------------------
  // Fit one side/charge category.
  // ------------------------------------------------------------

  bool FitCategory(
      TFile* inputFile,
      TFile* outputFile,
      const TString& category,
      const TString& pdfName,
      const bool firstPage,
      const bool lastPage)
  {
    TDirectory* inputDirectory =
        inputFile->GetDirectory(category);

    if (!inputDirectory)
    {
      std::cerr
          << "Directory not found: "
          << category
          << std::endl;

      return false;
    }

    TH2* sourceHistogram =
        GetInputHistogram(inputDirectory);

    if (!sourceHistogram)
    {
      std::cerr
          << "Could not find rDCA-vs-phi histogram in "
          << category
          << std::endl;

      return false;
    }

    outputFile->cd();

    TDirectory* outputDirectory =
        outputFile->GetDirectory(category);

    if (!outputDirectory)
    {
      outputDirectory =
          outputFile->mkdir(category);
    }

    outputDirectory->cd();

    TDirectory* sliceDirectory =
        outputDirectory->mkdir("slices");

    std::unique_ptr<TGraphErrors> graph(
        BuildSliceGraph(
            sourceHistogram,
            Form("g_rdca_slice_centers_%s", category.Data()),
            sliceDirectory));

    if (!graph)
    {
      std::cerr
          << "Could not construct slice graph for "
          << category
          << std::endl;

      return false;
    }

    const double xmin =
        -TMath::Pi();

    const double xmax =
        TMath::Pi();

    std::cout
        << "\n================================================\n"
        << "Category: " << category << "\n"
        << "Number of accepted slices: " << graph->GetN() << "\n"
        << "================================================"
        << std::endl;

    // ----------------------------------------------------------
    // Stage 1: A sin(x)
    // ----------------------------------------------------------

    TF1 fit1(
        Form("fit1_%s", category.Data()),
        "[0]*sin(x)",
        xmin,
        xmax);

    fit1.SetParName(0, "Main amplitude");
    fit1.SetParameter(0, 0.3);
    fit1.SetNpx(100000);

    graph->Fit(&fit1, "RQ0");

    // ----------------------------------------------------------
    // Stage 2: A sin(x) + offset
    // ----------------------------------------------------------

    TF1 fit2(
        Form("fit2_%s", category.Data()),
        "[0]*sin(x)+[1]",
        xmin,
        xmax);

    fit2.SetParNames(
        "Main amplitude",
        "Offset");

    fit2.SetParameter(
        0,
        fit1.GetParameter(0));

    fit2.SetParameter(
        1,
        0.0);

    fit2.SetNpx(100000);

    graph->Fit(&fit2, "RQ0");

    // ----------------------------------------------------------
    // Stage 3: A sin(x+phase) + offset
    // ----------------------------------------------------------

    TF1 fit3(
        Form("fit3_%s", category.Data()),
        "[0]*sin(x+[2])+[1]",
        xmin,
        xmax);

    fit3.SetParNames(
        "Main amplitude",
        "Offset",
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

    fit3.SetNpx(100000);

    for (int iteration = 0;
         iteration < 5;
         ++iteration)
    {
      graph->Fit(&fit3, "RQ0");
    }

    // ----------------------------------------------------------
    // Stage 4: common 12-fold amplitude, phase fixed initially.
    // ----------------------------------------------------------

    TF1 fit4(
        Form("fit4_%s", category.Data()),
        "[0]*sin(x+[2])+[1]+[3]*sin(12*x)",
        xmin,
        xmax);

    fit4.SetParNames(
        "Main amplitude",
        "Offset",
        "Main phase",
        "Common sector amplitude");

    fit4.SetParameter(
        0,
        fit3.GetParameter(0));

    fit4.SetParameter(
        1,
        fit3.GetParameter(1));

    fit4.SetParameter(
        2,
        fit3.GetParameter(2));

    fit4.SetParameter(
        3,
        0.1);

    fit4.SetParLimits(
        2,
        -TMath::Pi(),
        TMath::Pi());

    fit4.SetNpx(100000);

    for (int iteration = 0;
         iteration < NCommonFitIterations;
         ++iteration)
    {
      graph->Fit(&fit4, "RQ0");
    }

    // ----------------------------------------------------------
    // Stage 5: release common 12-fold phase.
    // ----------------------------------------------------------

    TF1 fit5(
        Form("fit5_%s", category.Data()),
        "[0]*sin(x+[2])+[1]+[3]*sin(12*x+[4])",
        xmin,
        xmax);

    fit5.SetParNames(
        "Main amplitude",
        "Offset",
        "Main phase",
        "Common sector amplitude",
        "Common sector phase");

    fit5.SetParameter(
        0,
        fit4.GetParameter(0));

    fit5.SetParameter(
        1,
        fit4.GetParameter(1));

    fit5.SetParameter(
        2,
        fit4.GetParameter(2));

    fit5.SetParameter(
        3,
        fit4.GetParameter(3));

    fit5.SetParameter(
        4,
        0.0);

    fit5.SetParLimits(
        2,
        -TMath::Pi(),
        TMath::Pi());

    fit5.SetParLimits(
        4,
        -TMath::Pi(),
        TMath::Pi());

    fit5.SetNpx(100000);

    for (int iteration = 0;
         iteration < NCommonFitIterations;
         ++iteration)
    {
      graph->Fit(&fit5, "RQ0");
    }

    // ----------------------------------------------------------
    // Final piecewise fit with 12 separate amplitudes.
    // ----------------------------------------------------------

    auto* finalFit =
        new TF1(
            Form("final_sector_fit_%s", category.Data()),
            FinalSectorFunction,
            xmin,
            xmax,
            NFinalParameters);

    finalFit->SetParName(0, "Main amplitude");
    finalFit->SetParName(1, "Offset");
    finalFit->SetParName(2, "Main phase");
    finalFit->SetParName(3, "Sector phase");
    finalFit->SetParName(4, "Second harmonic");

    finalFit->SetParameter(
        0,
        fit5.GetParameter(0));

    finalFit->SetParameter(
        1,
        fit5.GetParameter(1));

    finalFit->SetParameter(
        2,
        fit5.GetParameter(2));

    finalFit->SetParameter(
        3,
        fit5.GetParameter(4));

    finalFit->SetParameter(
        4,
        0.0);

    finalFit->SetParLimits(
        2,
        -TMath::Pi(),
        TMath::Pi());

    finalFit->SetParLimits(
        3,
        -TMath::Pi(),
        TMath::Pi());

    // Keep the second harmonic fixed initially.
    finalFit->FixParameter(4, 0.0);

    const double commonAmplitude =
        fit5.GetParameter(3);

    for (int sector = 0;
         sector < NSectors;
         ++sector)
    {
      const int parameter =
          5 + sector;

      finalFit->SetParName(
          parameter,
          Form("Sector %02d amplitude", sector));

      SetRelativeLimits(
          finalFit,
          parameter,
          commonAmplitude,
          RelativeStep);
    }

    finalFit->SetNpx(250000);

    for (int iteration = 0;
         iteration < NFinalFitIterations;
         ++iteration)
    {
      if (iteration > 0)
      {
        for (int sector = 0;
             sector < NSectors;
             ++sector)
        {
          const int parameter =
              5 + sector;

          const double previousValue =
              finalFit->GetParameter(parameter);

          SetRelativeLimits(
              finalFit,
              parameter,
              previousValue,
              RelativeStep);
        }
      }

      const TString options =
          (iteration == NFinalFitIterations - 1)
              ? "RSME"
              : "RQ0S";

      TFitResultPtr result =
          graph->Fit(
              finalFit,
              options);

      std::cout
          << "Final iteration "
          << iteration + 1
          << "/"
          << NFinalFitIterations
          << ", status = "
          << static_cast<int>(result)
          << ", chi2/NDF = "
          << finalFit->GetChisquare()
          << "/"
          << finalFit->GetNDF()
          << std::endl;
    }

    finalFit->SetLineColor(kRed + 1);
    finalFit->SetLineWidth(3);
    finalFit->SetNpx(250000);

    // ----------------------------------------------------------
    // Draw the slice-center graph.
    // ----------------------------------------------------------

    TCanvas canvas(
        Form("canvas_%s", category.Data()),
        category,
        1600,
        1000);

    canvas.SetLeftMargin(0.11);
    canvas.SetRightMargin(0.22);
    canvas.SetBottomMargin(0.12);
    canvas.SetTopMargin(0.10);

    graph->SetTitle(
        Form(
            "Slice-fit rDCA center vs atan2(p_{x},p_{y}) [%s];"
            "atan2(p_{x},p_{y});"
            "fitted rDCA center",
            category.Data()));

    graph->Draw("AP");

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

    outputDirectory->cd();

    graph->Write(
        graph->GetName(),
        TObject::kOverwrite);

    finalFit->Write(
        finalFit->GetName(),
        TObject::kOverwrite);

    canvas.Write(
        canvas.GetName(),
        TObject::kOverwrite);

    std::cout
        << "\nFinal parameters for "
        << category
        << "\nMain amplitude = "
        << finalFit->GetParameter(0)
        << " +/- "
        << finalFit->GetParError(0)
        << "\nOffset = "
        << finalFit->GetParameter(1)
        << " +/- "
        << finalFit->GetParError(1)
        << "\nMain phase = "
        << finalFit->GetParameter(2)
        << " +/- "
        << finalFit->GetParError(2)
        << "\nSector phase = "
        << finalFit->GetParameter(3)
        << " +/- "
        << finalFit->GetParError(3)
        << "\nChi2/NDF = "
        << finalFit->GetChisquare()
        << "/"
        << finalFit->GetNDF()
        << std::endl;

    for (int sector = 0;
         sector < NSectors;
         ++sector)
    {
      const int parameter =
          5 + sector;

      std::cout
          << "Sector "
          << sector
          << " amplitude = "
          << finalFit->GetParameter(parameter)
          << " +/- "
          << finalFit->GetParError(parameter)
          << std::endl;
    }

    outputFile->cd();

    return true;
  }
}

// ------------------------------------------------------------
// Main macro
// ------------------------------------------------------------

void FitRDCASlices(
    const char* inputFileName = "output/testpp_pcaz10.root",
    const char* outputPrefix = "outputfit/rdca_slice_fit")
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);

  gSystem->mkdir("output", true);

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
        << "Could not create output file: "
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
      << "\n================================================\n"
      << "Finished "
      << successfulFits
      << " of "
      << categories.size()
      << " categories.\n"
      << "Output ROOT file: "
      << outputRootName
      << "\n"
      << "Output PDF:       "
      << outputPdfName
      << "\n"
      << "================================================"
      << std::endl;
}