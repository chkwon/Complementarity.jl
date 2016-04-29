$Title Transportation model as equilibrium problem (TRANSMCP,SEQ=126)
$Ontext
   Dantzig's original transportation model (TRNSPORT) is
   reformulated as a linear complementarity problem.  We first
   solve the model with fixed demand and supply quantities, and
   then we incorporate price-responsiveness on both sides of the
   market.


Dantzig, G B, Chapter 3.3. In Linear Programming and Extensions.
Princeton University Press, Princeton, New Jersey, 1963.

[Source] http://www.gams.com/modlib/libhtml/transmcp.htm
$Offtext

  Sets
       i   canning plants   / seattle, san-diego /
       j   markets          / new-york, chicago, topeka / ;

  Parameters

       a(i)  capacity of plant i in cases (when prices are unity)
         /    seattle     350
              san-diego   600  /,

       b(j)  demand at market j in cases (when prices equal unity)
         /    new-york    325
              chicago     300
              topeka      275  /,

        esub(j)  price elasticity of demand (at prices equal to unity)
         /    new-york    1.5
              chicago     1.2
              topeka      2.0  /;

  Table d(i,j)  distance in thousands of miles
                    new-york       chicago      topeka
      seattle          2.5           1.7          1.8
      san-diego        2.5           1.8          1.4  ;

  Scalar f  freight in dollars per case per thousand miles  /90/ ;

  Parameter c(i,j)  transport cost in thousands of dollars per case ;

            c(i,j) = f * d(i,j) / 1000 ;

  Parameter pbar(j) reference price at demand node j;

  Positive variables
       w(i)             shadow price at supply node i,
       p(j)             shadow price at demand node j,
       x(i,j)           shipment quantities in cases;

  Equations
       supply(i)        supply limit at plant i,
       fxdemand(j)      fixed demand at market j,
       prdemand(j)      price-responsive demand at market j,
       profit(i,j)      zero profit conditions;

profit(i,j)..   w(i) + c(i,j)   =g= p(j);

supply(i)..     a(i) =g= sum(j, x(i,j));

fxdemand(j)..   sum(i, x(i,j))  =g=  b(j);

prdemand(j)..   sum(i, x(i,j))  =g=  b(j) * (pbar(j)/p(j))**esub(j);

*       declare models including specification of equation-variable
*       association:

  Model fixedqty / profit.x, supply.w, fxdemand.p/ ;
  Model equilqty / profit.x, supply.w, prdemand.p/ ;

*       initial estimate:

  p.l(j) = 1;
  w.l(i) = 1;

  Parameter report(*,*,*)  summary report;

  Solve fixedqty using mcp;
  report(i,j,"fixed") = x.l(i,j);
  report("price",j,"fixed") = p.l(j);
  report(i,"price","fixed") = w.l(i);

*       calibrate the demand functions:
  pbar(j) = p.l(j);

*       replicate the fixed demand equilibrium using flexible demand func:
  Solve equilqty using mcp;
  report(i,j,"flex") = x.l(i,j);
  report("price",j,"flex") = p.l(j);
  report(i,"price","flex") = w.l(i);

*       compute a counter-factual equilibrium using fixed demand func:
  c("seattle","chicago") = 0.5 * c("seattle","chicago");
  Solve fixedqty using mcp;
  report(i,j,"fixed CF") = x.l(i,j);
  report("price",j,"fixed CF") = p.l(j);
  report(i,"price","fixed CF") = w.l(i);

*       compute a counter-factual equilibrium using flexible demand func:
  Solve equilqty using mcp;
  report(i,j,"flex CF") = x.l(i,j);
  report("price",j,"flex CF") = p.l(j);
  report(i,"price","flex CF") = w.l(i);


  Display report;
execute_unload 'mcpReport.gdx', report;
