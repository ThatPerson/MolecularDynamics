import numpy as np
import sys

try:
	fn = sys.argv[1]
	out_fn = sys.argv[2]
except IndexError:
	print("Please pass filename")
	exit()

def movemol(buff, ofn):
	if (len(buff) == 0):
		return -1
	print(buff)
	num_atoms = int(buff[0])
	arr = np.zeros((3, num_atoms))
	current = 0
	for i in buff[2:]:
		k = i.strip().split()
		p = [float(l) for l in k[1:]]
		arr[:, current] = p
		current = current + 1
	com = np.mean(arr, axis=1)
	for i in buff:
		k = i.strip().split()
		if (len(k) == 1):
			ofn.write(i)
			continue
		p = np.array([float(l) for l in k[1:]])
		p = p - com
		ofn.write("%s\t%f\t%f\t%f\n" % (k[0], p[0], p[1], p[2]))
		
		

def is_int(k):
	try:
		int(k)
		return True
	except ValueError:
		return False

buff = []
with open(fn, "r") as f, open(out_fn, "w") as ofn:
	for l in f:
		k = l.split()
		if (len(k) == 1 and is_int(k[0].strip())):
			print("New")
			movemol(buff, ofn)
			buff = []
			buff.append(l)
			continue
		buff.append(l)
