#include "dogdefs.h"
#include <cmath>
#include "tensors.h"
#include "constants.h"
#include "Quadrature.h"

// ------------------------------------------------------------------------- //
// Set the 1D Gaussian quadrature weights and points.  These rules are designed
// to integrate over the interval [-1,1].
//
// This N-point quadrature rule
// integrates polynomials of degree 2*N-3 exactly.
//
// Returns:
// --------
//
//      w1d( 1:numpts )    - 1D quadrature weights
//      x1d( 1:numpts )    - 1D quadrature points
//
// Integration (over the canonical interval [-1,1]) is then given by
// 
//          \int_{-1}^1 f(x) dx \approx \sum_{i=1}^n f( x1d(i) ) * w1d( i ).
//
// See also: setGaussLobattoPoints1d
// ------------------------------------------------------------------------- //
void setGaussPoints1d(dTensor1& w1d, dTensor1& x1d)
{

    const int md2 = x1d.getsize()/2;   // used for very high-order quadrature

    // assert_eq( x1d.getsize(), w1d.getsize() );
    switch( x1d.getsize() )
    {
        case 0:
            break;

        case 1:
            w1d.set(1, 2.0 );

            x1d.set(1, 0.0 );
            break;

        case 2:
            w1d.set(1,  1.0 );
            w1d.set(2,  1.0 );

            x1d.set(1, -1.0/sq3 );
            x1d.set(2,  1.0/sq3 );
            break;

        case 3:
            w1d.set(1,  5.0/9.0 );
            w1d.set(2,  8.0/9.0 );
            w1d.set(3,  5.0/9.0 );

            x1d.set(1,  sq3/sq5 );
            x1d.set(2,  0.0 );
            x1d.set(3, -sq3/sq5 );
            break;

        case 4:
            w1d.set(1, (18.0 - sq3*sq10)/36.0 );
            w1d.set(2, (18.0 + sq3*sq10)/36.0 );
            w1d.set(3, w1d.get(2) );
            w1d.set(4, w1d.get(1) );

            x1d.set(1,  sqrt(3.0 + sqrt(4.8))/(-sq7) );
            x1d.set(2,  sqrt(3.0 - sqrt(4.8))/(-sq7) );
            x1d.set(3, -x1d.get(2) );
            x1d.set(4, -x1d.get(1) );       
            break;

        case 5:      
            w1d.set(1, (322.0 - 13.0*sq7*sq10)/900.0 );
            w1d.set(2, (322.0 + 13.0*sq7*sq10)/900.0 );
            w1d.set(3, 128.0/225.0 );
            w1d.set(4, w1d.get(2) );
            w1d.set(5, w1d.get(1) );

            x1d.set(1,  sqrt(5.0 + 2.0*sq10/sq7)/(-3.0) );
            x1d.set(2,  sqrt(5.0 - 2.0*sq10/sq7)/(-3.0) );
            x1d.set(3,  0.0 );
            x1d.set(4, -x1d.get(2) );
            x1d.set(5, -x1d.get(1) );
            break;

        case 6:

            w1d.set(1,  0.171324492379170);
            w1d.set(2,  0.360761573048139);
            w1d.set(3,  0.467913934572691);
            w1d.set(4,  0.467913934572691);
            w1d.set(5,  0.360761573048139);
            w1d.set(6,  0.171324492379170);


            x1d.set(1,  0.932469514203152);
            x1d.set(2,  0.661209386466265);
            x1d.set(3,  0.238619186083197);
            x1d.set(4, -x1d.get(3) );
            x1d.set(5, -x1d.get(2) );
            x1d.set(6, -x1d.get(1) );

            break;

        case 8:
            w1d.set(1, 0.1012285362903762591525314 );
            w1d.set(2, 0.2223810344533744705443560 );
            w1d.set(3, 0.3137066458778872873379622 );
            w1d.set(4, 0.3626837833783619829651504 );

            x1d.set(1, 0.9602898564975362316835609 );
            x1d.set(2, 0.7966664774136267395915539 );
            x1d.set(3, 0.5255324099163289858177390 );
            x1d.set(4, 0.1834346424956498049394761 );

            for (int k=1; k<=md2; k++)
            { w1d.set(md2+k,  w1d.get(md2+1-k) ); }      
            for (int k=1; k<=md2; k++)
            { x1d.set(md2+k, -x1d.get(md2+1-k) ); }
            break;

        case 10:
            w1d.set(1, 0.0666713443086881375935688 );
            w1d.set(2, 0.1494513491505805931457763 );
            w1d.set(3, 0.2190863625159820439955349 );
            w1d.set(4, 0.2692667193099963550912269 );
            w1d.set(5, 0.2955242247147528701738930 );

            x1d.set(1, 0.9739065285171717200779640 );
            x1d.set(2, 0.8650633666889845107320967 );
            x1d.set(3, 0.6794095682990244062343274 );
            x1d.set(4, 0.4333953941292471907992659 );
            x1d.set(5, 0.1488743389816312108848260 );

            for (int k=1; k<=md2; k++)
            { w1d.set(md2+k,  w1d.get(md2+1-k) ); }      
            for (int k=1; k<=md2; k++)
            { x1d.set(md2+k, -x1d.get(md2+1-k) ); }
            break;

        case 12:
            w1d.set(1, 0.0471753363865118271946160 );
            w1d.set(2, 0.1069393259953184309602547 );
            w1d.set(3, 0.1600783285433462263346525 );
            w1d.set(4, 0.2031674267230659217490645 );
            w1d.set(5, 0.2334925365383548087608499 );
            w1d.set(6, 0.2491470458134027850005624 );

            x1d.set(1, 0.9815606342467192506905491 );
            x1d.set(2, 0.9041172563704748566784659 );
            x1d.set(3, 0.7699026741943046870368938 );
            x1d.set(4, 0.5873179542866174472967024 );
            x1d.set(5, 0.3678314989981801937526915 );
            x1d.set(6, 0.1252334085114689154724414 );

            for (int k=1; k<=md2; k++)
            { w1d.set(md2+k,  w1d.get(md2+1-k) ); }      
            for (int k=1; k<=md2; k++)
            { x1d.set(md2+k, -x1d.get(md2+1-k) ); }
            break;

        case 14:
            w1d.set(1, 0.0351194603317518630318329 );
            w1d.set(2, 0.0801580871597602098056333 );
            w1d.set(3, 0.1215185706879031846894148 );
            w1d.set(4, 0.1572031671581935345696019 );
            w1d.set(5, 0.1855383974779378137417166 );
            w1d.set(6, 0.2051984637212956039659241 );
            w1d.set(7, 0.2152638534631577901958764 );

            x1d.set(1, 0.9862838086968123388415973 );
            x1d.set(2, 0.9284348836635735173363911 );
            x1d.set(3, 0.8272013150697649931897947 );
            x1d.set(4, 0.6872929048116854701480198 );
            x1d.set(5, 0.5152486363581540919652907 );
            x1d.set(6, 0.3191123689278897604356718 );
            x1d.set(7, 0.1080549487073436620662447 );

            for (int k=1; k<=md2; k++)
            { w1d.set(md2+k,  w1d.get(md2+1-k) ); }      
            for (int k=1; k<=md2; k++)
            { x1d.set(md2+k, -x1d.get(md2+1-k) ); }
            break;

        case 16:
            w1d.set(1, 0.0271524594117540948517806 );
            w1d.set(2, 0.0622535239386478928628438 );
            w1d.set(3, 0.0951585116824927848099251 );
            w1d.set(4, 0.1246289712555338720524763 );
            w1d.set(5, 0.1495959888165767320815017 );
            w1d.set(6, 0.1691565193950025381893121 );
            w1d.set(7, 0.1826034150449235888667637 );
            w1d.set(8, 0.1894506104550684962853967 );

            x1d.set(1, 0.9894009349916499325961542 );
            x1d.set(2, 0.9445750230732325760779884 );
            x1d.set(3, 0.8656312023878317438804679 );
            x1d.set(4, 0.7554044083550030338951012 );
            x1d.set(5, 0.6178762444026437484466718 );
            x1d.set(6, 0.4580167776572273863424194 );
            x1d.set(7, 0.2816035507792589132304605 );
            x1d.set(8, 0.0950125098376374401853193 );

            for (int k=1; k<=md2; k++)
            { w1d.set(md2+k,  w1d.get(md2+1-k) ); }      
            for (int k=1; k<=md2; k++)
            { x1d.set(md2+k, -x1d.get(md2+1-k) ); }
            break;

        case 18:
            w1d.set(1, 0.0216160135264833103133427 );
            w1d.set(2, 0.0497145488949697964533349 );
            w1d.set(3, 0.0764257302548890565291297 );
            w1d.set(4, 0.1009420441062871655628140 );
            w1d.set(5, 0.1225552067114784601845191 );
            w1d.set(6, 0.1406429146706506512047313 );
            w1d.set(7, 0.1546846751262652449254180 );
            w1d.set(8, 0.1642764837458327229860538 );
            w1d.set(9, 0.1691423829631435918406565 );

            x1d.set(1, 0.9915651684209309467300160 );
            x1d.set(2, 0.9558239495713977551811959 );
            x1d.set(3, 0.8926024664975557392060606 );
            x1d.set(4, 0.8037049589725231156824175 );
            x1d.set(5, 0.6916870430603532078748911 );
            x1d.set(6, 0.5597708310739475346078715 );
            x1d.set(7, 0.4117511614628426460359318 );
            x1d.set(8, 0.2518862256915055095889729 );
            x1d.set(9, 0.0847750130417353012422619 );

            for (int k=1; k<=md2; k++)
            { w1d.set(md2+k,  w1d.get(md2+1-k) ); }      
            for (int k=1; k<=md2; k++)
            { x1d.set(md2+k, -x1d.get(md2+1-k) ); }
            break;

        case 20:
            w1d.set(1,  0.0176140071391521183118620 );
            w1d.set(2,  0.0406014298003869413310400 );
            w1d.set(3,  0.0626720483341090635695065 );
            w1d.set(4,  0.0832767415767047487247581 );
            w1d.set(5,  0.1019301198172404350367501 );
            w1d.set(6,  0.1181945319615184173123774 );
            w1d.set(7,  0.1316886384491766268984945 );
            w1d.set(8,  0.1420961093183820513292983 );
            w1d.set(9,  0.1491729864726037467878287 );
            w1d.set(10, 0.1527533871307258506980843 );

            x1d.set(1,  0.9931285991850949247861224 );
            x1d.set(2,  0.9639719272779137912676661 );
            x1d.set(3,  0.9122344282513259058677524 );
            x1d.set(4,  0.8391169718222188233945291 );
            x1d.set(5,  0.7463319064601507926143051 );
            x1d.set(6,  0.6360536807265150254528367 );
            x1d.set(7,  0.5108670019508270980043641 );
            x1d.set(8,  0.3737060887154195606725482 );
            x1d.set(9,  0.2277858511416450780804962 );
            x1d.set(10, 0.0765265211334973337546404 );

            for (int k=1; k<=md2; k++)
            { w1d.set(md2+k,  w1d.get(md2+1-k) ); }      
            for (int k=1; k<=md2; k++)
            { x1d.set(md2+k, -x1d.get(md2+1-k) ); }
            break;

        case 50:
            x1d.set(1,  9.988664044200710e-01);
            x1d.set(2,  9.940319694320907e-01);
            x1d.set(3,  9.853540840480060e-01);
            x1d.set(4,  9.728643851066920e-01);
            x1d.set(5,  9.566109552428079e-01);
            x1d.set(6,  9.366566189448780e-01);
            x1d.set(7,  9.130785566557917e-01);
            x1d.set(8,  8.859679795236131e-01);
            x1d.set(9,  8.554297694299460e-01);
            x1d.set(10, 8.215820708593360e-01);
            x1d.set(11, 7.845558329003994e-01);
            x1d.set(12, 7.444943022260686e-01);
            x1d.set(13, 7.015524687068222e-01);
            x1d.set(14, 6.558964656854394e-01);
            x1d.set(15, 6.077029271849503e-01);
            x1d.set(16, 5.571583045146502e-01);
            x1d.set(17, 5.044581449074643e-01);
            x1d.set(18, 4.498063349740388e-01);
            x1d.set(19, 3.934143118975651e-01);
            x1d.set(20, 3.355002454194373e-01);
            x1d.set(21, 2.762881937795321e-01);
            x1d.set(22, 2.160072368760417e-01);
            x1d.set(23, 1.548905899981459e-01);
            x1d.set(24, 9.317470156008612e-02);
            x1d.set(25, 3.109833832718883e-02);
            for( int n=1; n <= 25; n++ ) 
            { x1d.set(25+n, -x1d.get(26-n) ); }

            w1d.set(1,  2.908622553155257e-03);
            w1d.set(2,  6.759799195745480e-03);
            w1d.set(3,  1.059054838365105e-02);
            w1d.set(4,  1.438082276148557e-02);
            w1d.set(5,  1.811556071348930e-02);
            w1d.set(6,  2.178024317012479e-02); 
            w1d.set(7,  2.536067357001242e-02); 
            w1d.set(8,  2.884299358053520e-02); 
            w1d.set(9,  3.221372822357796e-02); 
            w1d.set(10, 3.545983561514608e-02); 
            w1d.set(11, 3.856875661258761e-02); 
            w1d.set(12, 4.152846309014772e-02); 
            w1d.set(13, 4.432750433880325e-02); 
            w1d.set(14, 4.695505130394845e-02); 
            w1d.set(15, 4.940093844946636e-02); 
            w1d.set(16, 5.165570306958114e-02); 
            w1d.set(17, 5.371062188899628e-02); 
            w1d.set(18, 5.555774480621253e-02); 
            w1d.set(19, 5.718992564772843e-02); 
            w1d.set(20, 5.860084981322242e-02); 
            w1d.set(21, 5.978505870426547e-02); 
            w1d.set(22, 6.073797084177022e-02); 
            w1d.set(23, 6.145589959031677e-02); 
            w1d.set(24, 6.193606742068321e-02); 
            w1d.set(25, 6.217661665534725e-02); 
            for( int n=1; n <= 25; n++ ) 
            { w1d.set(25+n, w1d.get(26-n) ); }

            break;
        default:
            unsupported_value_error(x1d.getsize());
    }
}

// Same as above, but without the weights
void setGaussPoints1d(dTensor1& x1d)
{

    switch( x1d.getsize() )
    {
        case 0:
            break;

        case 1:
            x1d.set(1, 0.0 );
            break;

        case 2:
            x1d.set(1, -1.0/sq3 );
            x1d.set(2,  1.0/sq3 );
            break;

        case 3:
            x1d.set(1,  sq3/sq5 );
            x1d.set(2,  0.0 );
            x1d.set(3, -sq3/sq5 );
            break;

        case 4:
            x1d.set(1,  sqrt(3.0 + sqrt(4.8))/(-sq7) );
            x1d.set(2,  sqrt(3.0 - sqrt(4.8))/(-sq7) );
            x1d.set(3, -x1d.get(2) );
            x1d.set(4, -x1d.get(1) );       
            break;

        case 5:      
            x1d.set(1,  sqrt(5.0 + 2.0*sq10/sq7)/(-3.0) );
            x1d.set(2,  sqrt(5.0 - 2.0*sq10/sq7)/(-3.0) );
            x1d.set(3,  0.0 );
            x1d.set(4, -x1d.get(2) );
            x1d.set(5, -x1d.get(1) );
            break;

        case 6:
            x1d.set(1,  0.932469514203152);
            x1d.set(2,  0.661209386466265);
            x1d.set(3,  0.238619186083197);
            x1d.set(4, -x1d.get(3) );
            x1d.set(5, -x1d.get(2) );
            x1d.set(6, -x1d.get(1) );

            break;


        case 20:
            x1d.set(1,  0.993128599185095);
            x1d.set(2,  0.963971927277914);
            x1d.set(3,  0.912234428251326);
            x1d.set(4,  0.839116971822219);
            x1d.set(5,  0.746331906460151);
            x1d.set(6,  0.636053680726515);
            x1d.set(7,  0.510867001950827);
            x1d.set(8,  0.373706088715420);
            x1d.set(9,  0.227785851141645);
            x1d.set(10, 0.076526521133497);
            for( int n=1; n <= 10; n++ ) 
            { x1d.set(10+n, -x1d.get(11-n) ); }

            break;

        case 50:
            x1d.set(1,  9.988664044200710e-01);
            x1d.set(2,  9.940319694320907e-01);
            x1d.set(3,  9.853540840480060e-01);
            x1d.set(4,  9.728643851066920e-01);
            x1d.set(5,  9.566109552428079e-01);
            x1d.set(6,  9.366566189448780e-01);
            x1d.set(7,  9.130785566557917e-01);
            x1d.set(8,  8.859679795236131e-01);
            x1d.set(9,  8.554297694299460e-01);
            x1d.set(10, 8.215820708593360e-01);
            x1d.set(11, 7.845558329003994e-01);
            x1d.set(12, 7.444943022260686e-01);
            x1d.set(13, 7.015524687068222e-01);
            x1d.set(14, 6.558964656854394e-01);
            x1d.set(15, 6.077029271849503e-01);
            x1d.set(16, 5.571583045146502e-01);
            x1d.set(17, 5.044581449074643e-01);
            x1d.set(18, 4.498063349740388e-01);
            x1d.set(19, 3.934143118975651e-01);
            x1d.set(20, 3.355002454194373e-01);
            x1d.set(21, 2.762881937795321e-01);
            x1d.set(22, 2.160072368760417e-01);
            x1d.set(23, 1.548905899981459e-01);
            x1d.set(24, 9.317470156008612e-02);
            x1d.set(25, 3.109833832718883e-02);
            for( int n=1; n <= 25; n++ ) 
            { x1d.set(25+n, -x1d.get(26-n) ); }

            break;


        default:
            unsupported_value_error(x1d.getsize());
    }
}

// ------------------------------------------------------------------------- //
// Set the 1D Gauss-Lobatto quadrature points.  These rules are designed
// to integrate over the interval [-1,1].  
// 
// This N-point quadrature rule
// integrates polynomials of degree 2*N-3 exactly.
//
// Returns:
// --------
//
//      x1d( 1:numpts )    - 1D quadrature points
//      w1d( 1:numpts )    - 1D quadrature weights
//
// See also: setGaussWeightsAndPoints
// ------------------------------------------------------------------------- //
void setGaussLobattoPoints1d( dTensor1& w1d, dTensor1& x1d)
{

    switch( x1d.getsize() )
    {
        case 2:
            x1d.set(1, -1.0 );
            x1d.set(2,  1.0 );

            w1d.set(1, 1.0 );
            w1d.set(2, 1.0 );
            break;

        case 3:
            x1d.set(1,-1.0);
            x1d.set(2, 0.0);
            x1d.set(3, 1.0);

            w1d.set(1, 1.0/3.0 );
            w1d.set(2, 4.0/3.0 );
            w1d.set(3, 1.0/3.0 );
            break;

        case 4:
            x1d.set(1, -1.0 );
            x1d.set(2, -0.447213595499958 );
            x1d.set(3,  0.447213595499958 );
            x1d.set(4, 1.0 );


            w1d.set(1, 0.166666666666667 );
            w1d.set(2, 0.833333333333333 );
            w1d.set(3, w1d.get(2) );
            w1d.set(4, w1d.get(1) );
            break;

        case 5:
            x1d.set(1, -1.0 );
            x1d.set(2, -0.654653670707977 );
            x1d.set(3,  0.0 );
            x1d.set(4, -x1d.get(2) );
            x1d.set(5, -x1d.get(1) );

            w1d.set(1, 0.10 );
            w1d.set(2, 0.544444444444444 );
            w1d.set(3, 0.711111111111111 );
            w1d.set(4, w1d.get(2) );
            w1d.set(5, w1d.get(1) );
            break;

        case 6:
            x1d.set(1, -1.0 );
            x1d.set(2, -0.765055323929465 );
            x1d.set(3, -0.285231516480645 );
            x1d.set(4, -x1d.get(3) );
            x1d.set(5, -x1d.get(2) );
            x1d.set(6, -x1d.get(1) );

            w1d.set(1, 0.066666666666667 );
            w1d.set(2, 0.378474956297847 ); 
            w1d.set(3, 0.554858377035486 );
            w1d.set(4, w1d.get(3) );
            w1d.set(5, w1d.get(2) );
            w1d.set(6, w1d.get(1) );
            break;

        case 7:      
            x1d.set(1,  -1);
            x1d.set(2,  -8.302238962785670e-01);
            x1d.set(3,  -4.688487934707142e-01);
            x1d.set(4,   0);
            x1d.set(5, -x1d.get(3) );
            x1d.set(6, -x1d.get(2) );
            x1d.set(7, -x1d.get(1) );

            w1d.set(1, 0.047619047619048 );
            w1d.set(2, 0.276826047361566 );
            w1d.set(3, 0.431745381209863 );
            w1d.set(4, 0.487619047619048 );
            w1d.set(5, w1d.get(3) );
            w1d.set(6, w1d.get(2) );
            w1d.set(7, w1d.get(1) );

            break;

        default:
            exit(1);
    }
}


void setGaussLobattoPoints1d(dTensor1& x1d)
{
    switch(x1d.getsize())
    {
        case 2:
            x1d.set(1, -1.0 );
            x1d.set(2,  1.0 );
            break;

        case 3:
            x1d.set(1,-1.0);
            x1d.set(2, 0.0);
            x1d.set(3, 1.0);
            break;

        case 4:
            x1d.set(1, -1.0 );
            x1d.set(2, -0.447213595499958 );
            x1d.set(3,  0.447213595499958 );
            x1d.set(4, 1.0 );
            break;

        case 5:
            x1d.set(1, -1.0 );
            x1d.set(2, -0.654653670707977 );
            x1d.set(3,  0.0 );
            x1d.set(4, -x1d.get(2) );
            x1d.set(5, -x1d.get(1) );
            break;

        case 6:
            x1d.set(1, -1.0 );
            x1d.set(2, -0.765055323929465 );
            x1d.set(3, -0.285231516480645 );
            x1d.set(4, -x1d.get(3) );
            x1d.set(5, -x1d.get(2) );
            x1d.set(6, -x1d.get(1) );
            break;

        case 7:      
            x1d.set(1,  -1);
            x1d.set(2,  -8.302238962785670e-01);
            x1d.set(3,  -4.688487934707142e-01);
            x1d.set(4,   0);
            x1d.set(5, -x1d.get(3) );
            x1d.set(6, -x1d.get(2) );
            x1d.set(7, -x1d.get(1) );
            break;

        default:
            unsupported_value_error(x1d.getsize());
    }
}


