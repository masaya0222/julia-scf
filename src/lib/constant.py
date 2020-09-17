import os
lib = os.path.dirname(__file__) + "/autogen"
am = 2
new_am = am*2 #10
old_am = 0
opt_am = new_am
am_letter = ["0", "p", "d", "f", "g", "h", "i", "k", "l", "m","n","o","q","r","t","u","v","w","x","y","z"]
number = ["zero", "one", "two", "three", "four", "five","six","seven","eight","nine","ten","eleven",
			"twelve","thirteen","fourteen","fifteen","sixteen","seventeen","eighteen","nineteen","twenty"]

MAX_AM = 22

def io(i):
    return (i*(i+1))//2

def hash(a, b):
    c = [0,0]
    if (b[0]):
        i = b[0]-a[0][0]
        c[0] = i + io(i) -a[0][1]
    if b[1]:
        i = b[1]-a[1][0]
        c[1] = i+io(i)-a[1][1]
    return c[0]*io(b[1]+1)+c[1]