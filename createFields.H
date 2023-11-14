#include "createRDeltaT.H"

Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"

Info<< "Constructing twoPhaseMixtureThermo\n" << endl;
twoPhaseMixtureThermo mixture(U, phi);

volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());

Info<< "Constructing compressiblePhaseChangeMixture\n" <<endl;

autoPtr<compressiblePhaseChangeMixture> fluidPtr =
    compressiblePhaseChangeMixture::New(mixture, mesh);
compressiblePhaseChangeMixture & fluid = fluidPtr();

Info<< "Reading thermophysical properties\n" << endl;









const volScalarField& rho1 = fluid.rho1();
const volScalarField& rho2 = fluid.rho2();

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    alpha1*rho1 + alpha2*rho2
);

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);


/*volScalarField vaporMassFraction
(
    IOobject
    (
        "vaporMassFraction",    // 组分名称
        runTime.timeName(),     // 当前时间
        mesh,                   // 网格
        IOobject::MUST_READ,      // 不从文件读取
        IOobject::AUTO_WRITE    // 自动写入
    ),
    mesh,
    dimensionedScalar("vaporMassFraction", dimless, 0.0)  // 初始化为零
);*/
volScalarField vaporMassFraction
(
    IOobject
    (
        "vaporMassFraction",    // 组分名称
        runTime.timeName(),     // 当前时间
        mesh,                   // 网格
        IOobject::MUST_READ,    // 从文件读取
        IOobject::AUTO_WRITE    // 自动写入
    ),
    mesh
);


volScalarField Su1
(
    IOobject
    (
        "Su1",   
        runTime.timeName(),     // 当前时间
        mesh,                   // 网格
        IOobject::NO_READ,      // 不从文件读取
        IOobject::AUTO_WRITE    // 自动写入
    ),
    mesh,
    dimensionedScalar("Su1", dimensionSet(0, 0, -1, 0, 0, 0, 0), 0.0)  // 初始化为零
);
volScalarField Su2
(
    IOobject
    (
        "Su2",   
        runTime.timeName(),     // 当前时间
        mesh,                   // 网格
        IOobject::NO_READ,      // 不从文件读取
        IOobject::AUTO_WRITE    // 自动写入
    ),
    mesh,
    dimensionedScalar("Su2", dimensionSet(0, 0, -1, 0, 0, 0, 0), 0.0)  // 初始化为零
);

scalar D = readScalar(transportProperties.lookup("D"));


dimensionedScalar myD
(
    "myD",
    dimensionSet(0, 2, -1, 0, 0, 0, 0),
    D
);

dimensionedScalar pMin
(
    "pMin",
    dimPressure,
    mixture
);

mesh.setFluxRequired(p_rgh.name());
mesh.setFluxRequired(alpha1.name());


#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"


// Mass flux
// Initialisation does not matter because rhoPhi is reset after the
// alpha1 solution before it is used in the U equation.
surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rho)*phi
);

volScalarField dgdt(alpha1*fvc::div(phi));

#include "createAlphaFluxes.H"
Foam::isoAdvection advector(alpha1,phi,U);
// Construct compressible turbulence model
compressibleInterPhaseTransportModel turbulence
(
    rho,
    U,
    phi,
    rhoPhi,
    alphaPhi10,
    mixture
);

#include "createK.H"

#include "createMRF.H"
#include "createFvOptions.H"
