tmp<volScalarField> tdmdtNet = fluid.dmdtNet();
const volScalarField& dmdt12 = tdmdtNet();
tmp<volScalarField> tcoeffs(fluid.coeffs());
const volScalarField& coeffs = tcoeffs();

forAll(dmdt12, celli)
{

        scalar cdmdt = dmdt12[celli];
        scalar mDot = 1.0/rho1[celli]*cdmdt ; // 相变时的质量转移速率
         scalar ccoeffs = coeffs[celli];
         scalar alpha1Limited = max(min(alpha1[celli],1.0),0.0);
        if (mDot < 0) // 如果蒸发，则在气相添加质量分数
    {
        Su1[celli] += mDot*ccoeffs*alpha1Limited; // 质量分数的源项
       // Su2[celli] -= mDot*ccoeffs*alpha1Limited;
    }
}


// 蒸汽的输运方程
fvScalarMatrix Y1
(
    fvm::ddt(vaporMassFraction) + fvm::div(phi, vaporMassFraction)
  - fvm::laplacian((1-alpha1)*myD, vaporMassFraction) == Su1
);
Y1.relax();
Y1.solve();
vaporMassFraction.correctBoundaryConditions();

/*// 空气的输运方程
fvScalarMatrix Y2
(
    fvm::ddt(airMassFraction) + fvm::div(phi, airMassFraction)
  - fvm::laplacian(alpha1*myD, airMassFraction) == Su2
);
Y2.relax();
Y2.solve();
airMassFraction.correctBoundaryConditions();

forAll(vaporMassFraction,celli)
{

vaporMassFraction[celli] = max(min(vaporMassFraction[celli],1.0),0.0);
airMassFraction[celli] = 1-vaporMassFraction[celli];
}*/
