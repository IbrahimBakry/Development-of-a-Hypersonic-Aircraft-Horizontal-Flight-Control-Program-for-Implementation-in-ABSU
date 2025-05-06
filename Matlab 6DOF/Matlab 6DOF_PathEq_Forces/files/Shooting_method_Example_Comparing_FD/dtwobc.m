function res = dtwobc(ya,yb)
res = [ya(4) - 30000; 
       ya(5) - 33;
       ya(6) - 55;
       ya(7) - 80000;
       yb(4) - 30000;
       yb(5) - 45/57.3;
       yb(6) + 73/57.3];