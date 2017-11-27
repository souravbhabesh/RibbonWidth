f(x)= a*x**b
#f(x) = 1/(A2 * x**2 + A4 * x**4)
a = 10
b = -1
#A2=1
#A4=1
fit [3:25] f(x) "WIDTH_L101_W21.txt" u 1:2:3 via a,b
ti = sprintf("%.4fx**%.4f", a, b)
#set logscale x
#set logscale y
plot "WIDTH_L101_W21.txt" u 1:2:3 with errorlines lt rgb "#ff0000" title "L101 W21", \
f(x) with lines lt rgb "#ff00ff" title ti
