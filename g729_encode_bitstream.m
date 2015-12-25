function bits = g729_encode_bitstream(pkt)

  % format contained in Table (8)

  bits = '';
 
  % encode LSF params L0,L1,L2,L3
  bits(1:18)  = encode_lsf(pkt.lsfIdx);
  
  % encode 1st subframe  
  % encode pitch P1
  if( pkt.sf{1}.T_int <= 85)
    P1 = 3*(pkt.sf{1}.T_int - 19) + pkt.sf{1}.T_frac - 1;
  else
    P1 = pkt.sf{1}.T_int - 85 + 197;
  end
  
  P1_bits     = dec2bin(P1,8);
  bits(19:26) = P1_bits;
  
  % encode parity bit
  P0          = mod(sum(arrayfun(@bin2dec,P1_bits)),2);
  bits(27)    = dec2bin(P0);
  
  % encode fixed codebook first subframe, C1  
  bits(28:40) = encode_fcb(pkt.sf{1}.iPulse);
  
  % encode signs of fixed codebook first subframe, S1
  bits(41:44) = encode_signs(pkt.sf{1}.sPulse);
  
  % encode gain codebook stage 1 (GA1) and stage 2 (GB1)
  bits(45:51) = encode_gain(pkt.sf{1}.gainIdx);
   
  % encode second subframe 
  % encode pitch delay P2
  
  % reflect pitch offset, starting with [T1 - 5 2/3] as 0
  % end with [T1 + 4 1/3] as 31.
  
  % account for fact that starting range can't be less than 19 1/3
  Tstart = (3*pkt.sf{1}.T_int - 17);
  Tstart = max(Tstart,58);
  deltaT = (3*pkt.sf{2}.T_int + pkt.sf{2}.T_frac) - Tstart;
                                       
  assert(deltaT >= 0 && deltaT <= 31);  % validate fits in 5 bits
  bits(52:56) = dec2bin(deltaT,5);
  
  % encode fixed codebook second subframe, C2  
  bits(57:69) = encode_fcb(pkt.sf{2}.iPulse);

  % encode signs of fixed codebook second subframe, S2
  bits(70:73) = encode_signs(pkt.sf{2}.sPulse);

  % encode gain codebook stage 1 (GA2) and stage 2 (GB2) 
  bits(74:80) = encode_gain(pkt.sf{2}.gainIdx); 

end


function S = encode_signs(sPulse)

  % [s3 s2 s1 s0] in order
  sBits = arrayfun(@dec2bin,sPulse > 0);
  S = fliplr(sBits);

end

function G = encode_gain(gainIdx)

  G = '';

  G(1:3)   = dec2bin(gainIdx.iA,3);
  G(4:7)   = dec2bin(gainIdx.iB,4);   
end

function L = encode_lsf(lsfIdx)

  L = '';
  
  % encode LSF params L0,L1,L2,L3
  L(1)     = dec2bin(lsfIdx.i0,1);
  L(2:8)   = dec2bin(lsfIdx.i1,7);
  L(9:13)  = dec2bin(lsfIdx.i2,5);
  L(14:18) = dec2bin(lsfIdx.i3,5);  
end

function C = encode_fcb(iPulse)

  C = '';
  i           = floor(iPulse/5);              % convert pulse positions to indices 0:7
  C(1:3)      = dec2bin(i(3+(1)),3);
  C(4)        = dec2bin(mod(iPulse(3+(1)),5) == 4);  % check if m3 in track [3:5:38] or [4:5:39]
  C(5:7)      = dec2bin(i(2+(1)),3);
  C(8:10)     = dec2bin(i(1+(1)),3);
  C(11:13)    = dec2bin(i(0+(1)),3);
  
end