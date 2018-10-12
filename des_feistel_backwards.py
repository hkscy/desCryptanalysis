# Chris Hicks 08-10-2018
# Hypothesis: The 32-bit half-block input to the DES Feistel function can 
#             be uniquely determined if the output from the function and the
#             round subkey are both known.
#
# Test: Show that the half block can be uniquely determined by executing the 
#       DES Feistel function with the round-key known and fixed and the output
#		known and fixed
#
# Results: The input_block entropy can be reduced to <=8 (Seen 5 max) from 2^32
#
# Hypothesis 2: If there are more than one valid input block, they are collisions
#
# Test: Execute function on all candidates and evaluate output equality
#
# Results: All inbut blocks retrieved by first test are colliding

import des_header
import random
import itertools

# Return DES sbox inputs with row bits added
def add_row_bits(sbox_inputs):
	in0 = sbox_inputs[0]<<1 		# 00 row bits for input 0
	in1 = sbox_inputs[1]<<1|0x1		# 01 row bits for input 1
	in2 = sbox_inputs[2]<<1|0x20	# 10 row bits for input 2
	in3 = sbox_inputs[3]<<1|0x21	# 11 row bits for input 3
	return [in0,in1,in2,in3]

def main():
	inr = random.getrandbits(32)
	print('Random DES Feistel input: ' + str(inr))

	e = 0
	for bidx in range(48):
		b = (inr >> 31-des_header.E[bidx])&0x1
		e |= (b << 47-bidx)
	print('E = ' + str(e) + ' ' + str(len(bin(e)[2:])))
	#print([bin(e)[2:].zfill(48)[i:i+6] for i in range(0, len(bin(e)[2:].zfill(48)), 6)])

	# Apply S-Boxes to E
	sixbit_mask = [0xfc0000000000,0x3f000000000,0xfc0000000,0x3f000000,0xfc0000,0x3f000,0xfc0,0x3f]
	foursixbit_mask = [0x780000000000,0x1e000000000,0x780000000,0x1e000000,0x780000,0x1e000,0x780,0x1e]
	
	# First and last bits determine row in 2D sbox table
	row = ((e&sixbit_mask[7])&0x20)>>4 | (e&sixbit_mask[7]&0x1)
	s8 = des_header.S_BOX[7][row][(e&foursixbit_mask[7])>>1]

	row = (((e&sixbit_mask[6])>>6)&0x20)>>4 | (((e&sixbit_mask[6])>>6)&0x1)
	s7 = des_header.S_BOX[6][row][(e&foursixbit_mask[6])>>7]

	row = (((e&sixbit_mask[5])>>12)&0x20)>>4 | (((e&sixbit_mask[5])>>12)&0x1)
	s6 = des_header.S_BOX[5][row][(e&foursixbit_mask[5])>>13]

	row = (((e&sixbit_mask[4])>>18)&0x20)>>4 | (((e&sixbit_mask[4])>>18)&0x1)
	s5 = des_header.S_BOX[4][row][(e&foursixbit_mask[4])>>19]

	row = (((e&sixbit_mask[3])>>24)&0x20)>>4 | (((e&sixbit_mask[3])>>24)&0x1)
	s4 = des_header.S_BOX[3][row][(e&foursixbit_mask[3])>>25]

	row = (((e&sixbit_mask[2])>>30)&0x20)>>4 | (((e&sixbit_mask[2])>>30)&0x1)
	s3 = des_header.S_BOX[2][row][(e&foursixbit_mask[2])>>31]

	row = (((e&sixbit_mask[1])>>36)&0x20)>>4 | (((e&sixbit_mask[1])>>36)&0x1)
	s2 = des_header.S_BOX[1][row][(e&foursixbit_mask[1])>>37]

	row = (((e&sixbit_mask[0])>>42)&0x20)>>4 | (((e&sixbit_mask[0])>>42)&0x1)
	s1 = des_header.S_BOX[0][row][(e&foursixbit_mask[0])>>43]

	sbox_out = s1<<28|s2<<24|s3<<20|s4<<16|s5<<12|s6<<8|s7<<4|s8
	print('sbox output = ' + str(sbox_out))
	#print(bin(sbox_out)[2:].zfill(32))
	#print([bin(sbox_out)[2:].zfill(32)[i:i+4] for i in range(0, len(bin(sbox_out)[2:].zfill(32)), 4)])


	# Apply final permutation
	#p = des_header.apply_pbox_bitwise(sbox_out, des_header.P, 32)
	p = 0
	for bidx in range(32):
		b = (sbox_out >> 31-des_header.P[bidx])&0x1
		p |= (b << 31-bidx)

	print('Final output: ' + str(p))
	print()

	#Now try re-create input

	# Apply P inverse to determine output from sboxes S1-S8
	#p_inv = des_header.apply_pbox_bitwise(out, des_header.P_inv, 32)
	p_inv = 0
	for bidx in range(32):
		b = (p >> 31-des_header.P_inv[bidx])&0x1
		p_inv |= (b << 31-bidx)

	#print(bin(out)[2:].zfill(32))
	#print(bin(p_inv)[2:].zfill(32))

	print('P inverse: ' + str(p_inv))

	# Each 4 bits in P-inverse is output from SBoxes S1-S8, but each Sbox 
	# takes 6 bits as input so for each output value there are four possible
	# inputs which could have caused the output.

	# Mask for sbox output bits
	mask = [0xf0000000,0xf000000,0xf00000,0xf0000,0xf000,0xf00,0xf0,0xf]

	sbox_inputs = [0]*8
	# For each DES SBOX
	for sbox in range(7,-1,-1):
		# Get relevant four bits from p_inv
		sbox_out = (p_inv&mask[sbox])>>((7-sbox)*4)
		#print('                S' + str(sbox+1) + ' output: ' + str (bin(sbox_out)[2:].zfill(4)))
		
		# Retrieve candidate inputs
		sbox_inputs[sbox] = des_header.SBoxInv[sbox][sbox_out]
		#print(des_header.SBoxInv[sbox][sbox_out])
		sbox_inputs[sbox] = add_row_bits(sbox_inputs[sbox])
		#print([bin(k) for k in sbox_inputs[sbox]])

	# Uncomment to print sbox candidate input guesses
	#for si, sbox_input in enumerate(sbox_inputs):
	#	print('Candidate S' + str(si+1) + ' inputs: ' + str([bin(s)[2:].zfill(6) for s in sbox_input]))	

	# We want to reduce the total search space by looking for conflicts/impossibility results in the candidates
	# lazy method:
	candidates = []
	for s8 in sbox_inputs[7]:
		for s7 in sbox_inputs[6]:
			for s6 in sbox_inputs[5]:
				for s5 in sbox_inputs[4]:
					for s4 in sbox_inputs[3]:
						for s3 in sbox_inputs[2]:
							for s2 in sbox_inputs[1]:
								for s1 in sbox_inputs[0]:
									# Expansion permutation for each six bits is as follows
									# 32	1	2	3	4	5	# (32)s1.2==s8.6 and  (1)s1.1==s8.5
									# 4		5	6	7	8	9	#  (8)s2.2==s1.6 and  (4)s2.1==s1.5
									# 8		9	10	11	12	13  #  (9)s3.2==s2.6 and (12)s3.1==s2.5
									# 12	13	14	15	16	17  # (13)s4.2==s3.6 and (16)s4.1==s3.5
									# 16	17	18	19	20	21  # ...
									# 20	21	22	23	24	25  # etc
									# 24	25	26	27	28	29  #
									# 28	29	30	31	32	1	#
									# Because of expansion function certain properties constrain the candidates
									# input bit 1 goes to output bit 2 and 48: s1.2==s8.6
									if (s1&0x10)>>4==(s8&0x1) and (s1&0x20)>>5==(s8&0x2)>>1:
										if (s2&0x10)>>4==(s1&0x1) and (s2&0x20)>>5==(s1&0x2)>>1:
											if (s3&0x10)>>4==(s2&0x1) and (s3&0x20)>>5==(s2&0x2)>>1:
												if (s4&0x10)>>4==(s3&0x1) and (s4&0x20)>>5==(s3&0x2)>>1:
													if (s5&0x10)>>4==(s4&0x1) and (s5&0x20)>>5==(s4&0x2)>>1:
														if (s6&0x10)>>4==(s5&0x1) and (s6&0x20)>>5==(s5&0x2)>>1:
															if (s7&0x10)>>4==(s6&0x1) and (s7&0x20)>>5==(s6&0x2)>>1:
																if (s8&0x10)>>4==(s7&0x1) and (s8&0x20)>>5==(s7&0x2)>>1:
																	candidates += [(s1)<<42|(s2)<<36|(s3)<<30|s4<<24|(s5)<<18|(s6)<<12|(s7)<<6|(s8)]

	# We're done here -- we should have 16 or less candidates remaining

	print('Found ' + str(len(candidates)) + ' candidate inputs:\n')
	inputs = []
	for c in candidates:
		# Apply inverse expansion to E	
		e_inv = 0
		for bidx in range(32):
			b = (c >> 47-des_header.E_inv[bidx])&0x1
			e_inv |= (b << 31-bidx)

		inputs += [e_inv]
		if e_inv==inr:
			print('Holy smokes: ' + str(e_inv))
		else:
			print(e_inv)	

	for inBlock in inputs:
		print()
		e = 0
		for bidx in range(48):
			b = (inBlock >> 31-des_header.E[bidx])&0x1
			e |= (b << 47-bidx)
		print('E = ' + str(e) + ' ' + str(len(bin(e)[2:])))
		#print([bin(e)[2:].zfill(48)[i:i+6] for i in range(0, len(bin(e)[2:].zfill(48)), 6)])

		# Apply S-Boxes to E
		sixbit_mask = [0xfc0000000000,0x3f000000000,0xfc0000000,0x3f000000,0xfc0000,0x3f000,0xfc0,0x3f]
		foursixbit_mask = [0x780000000000,0x1e000000000,0x780000000,0x1e000000,0x780000,0x1e000,0x780,0x1e]
		
		# First and last bits determine row in 2D sbox table
		row = ((e&sixbit_mask[7])&0x20)>>4 | (e&sixbit_mask[7]&0x1)
		s8 = des_header.S_BOX[7][row][(e&foursixbit_mask[7])>>1]

		row = (((e&sixbit_mask[6])>>6)&0x20)>>4 | (((e&sixbit_mask[6])>>6)&0x1)
		s7 = des_header.S_BOX[6][row][(e&foursixbit_mask[6])>>7]

		row = (((e&sixbit_mask[5])>>12)&0x20)>>4 | (((e&sixbit_mask[5])>>12)&0x1)
		s6 = des_header.S_BOX[5][row][(e&foursixbit_mask[5])>>13]

		row = (((e&sixbit_mask[4])>>18)&0x20)>>4 | (((e&sixbit_mask[4])>>18)&0x1)
		s5 = des_header.S_BOX[4][row][(e&foursixbit_mask[4])>>19]

		row = (((e&sixbit_mask[3])>>24)&0x20)>>4 | (((e&sixbit_mask[3])>>24)&0x1)
		s4 = des_header.S_BOX[3][row][(e&foursixbit_mask[3])>>25]

		row = (((e&sixbit_mask[2])>>30)&0x20)>>4 | (((e&sixbit_mask[2])>>30)&0x1)
		s3 = des_header.S_BOX[2][row][(e&foursixbit_mask[2])>>31]

		row = (((e&sixbit_mask[1])>>36)&0x20)>>4 | (((e&sixbit_mask[1])>>36)&0x1)
		s2 = des_header.S_BOX[1][row][(e&foursixbit_mask[1])>>37]

		row = (((e&sixbit_mask[0])>>42)&0x20)>>4 | (((e&sixbit_mask[0])>>42)&0x1)
		s1 = des_header.S_BOX[0][row][(e&foursixbit_mask[0])>>43]

		sbox_out = s1<<28|s2<<24|s3<<20|s4<<16|s5<<12|s6<<8|s7<<4|s8
		print('sbox output = ' + str(sbox_out))
		#print(bin(sbox_out)[2:].zfill(32))
		#print([bin(sbox_out)[2:].zfill(32)[i:i+4] for i in range(0, len(bin(sbox_out)[2:].zfill(32)), 4)])


		# Apply final permutation
		#p = des_header.apply_pbox_bitwise(sbox_out, des_header.P, 32)
		p = 0
		for bidx in range(32):
			b = (sbox_out >> 31-des_header.P[bidx])&0x1
			p |= (b << 31-bidx)

		print('Final output: ' + str(p))
		print()

	return 0

if __name__ == '__main__':
	main()