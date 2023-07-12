      function fxmngau(x)
      common /pawpar/ par(2)
      fxmngau = par(1)*exp(-0.5*((x-18.903)/par(2))**2)
      end
