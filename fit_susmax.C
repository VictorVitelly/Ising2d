void fit_susmax() {
    TCanvas *c = new TCanvas("c", "Susceptibility Max Fit", 800, 600);
    c->SetGrid();

    // Read file: x = L, y = susmax
    TGraphErrors *gr = new
    TGraphErrors("data/tc.dat", "%lg %lg %lg");
    gr->SetTitle("Max Susceptibility vs L;L;#chi_{max}");
    gr->SetMarkerStyle(20);
    gr->SetMarkerSize(1.0);
    gr->SetMarkerColor(kBlue+1);
    gr->SetLineColor(kBlue+1);
    gr->Draw("AP");

    // Fit function: a*x^b + c
    TF1 *fitFunc = new TF1("fitFunc", "[0]",1,3);
    fitFunc->SetParameters(1);  // initial guesses
    fitFunc->SetLineColor(kRed);
    fitFunc->SetLineWidth(2);

    gr->Fit(fitFunc, "R");
    fitFunc->Draw("SAME");

    // Extract parameters and errors
    double a     = fitFunc->GetParameter(0);
    double aerr  = fitFunc->GetParError(0);

    double chi2 = fitFunc->GetChisquare();
    double ndf  = fitFunc->GetNDF();
    double chi2_ndf = chi2 / ndf;

    // Text box with results
    TPaveText *pt = new TPaveText(0.60, 0.65, 0.88, 0.88, "NDC");
    pt->SetFillColorAlpha(0, 0.4);
    pt->SetBorderSize(1);
    pt->SetTextSize(0.035);
    pt->AddText(Form("a = %.3g #pm %.3g", a, aerr));
    pt->AddText(Form("#chi^{2}/ndf = %.2f", chi2_ndf));
    pt->Draw();

    c->SaveAs("susmax_fit.png");

    printf("\nFit results:\n");
    printf("a = %g ± %g\n", a, aerr);
    printf("chi^2/ndf = %g\n", chi2_ndf);
}
