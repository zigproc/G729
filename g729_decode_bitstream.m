function pkt = g729_decode_bitstream(bits)

  % format contained in Table (8)

  pkt = [];
 
  % encode LSF params L0,L1,L2,L3
  pkt.lsfIdx  = decode_lsf(bits(1:18));
  
  % decode 1st subframe  
  % decode pitch P1
  P1_bits = bits(19:26);
  P1 = bin2dec(P1_bits);
  if( P1 > 197)
    T_int = P1 - 197 + 85;
    T_frac = 0;
  else
    % 58 represents 19 1/3 (using 1/3 resolution)
    tmp = P1 + 58;
    T_int = floor(tmp/3);
    T_frac = mod(tmp,3);
    if( T_frac == 2)
      T_int += 1;
      T_frac -= 3;
    end    
  end
  
  T3 = 3*T_int;
  
  pkt.sf{1}.T_int = T_int;
  pkt.sf{1}.T_frac = T_frac;
  
  % check parity bit
  P0 = bin2dec(bits(27));
  parityRx = mod(sum(arrayfun(@bin2dec,P1_bits)),2);
  assert( P0 == parityRx );
  
  % decode fixed codebook first subframe, C1  
  pkt.sf{1}.iPulse = decode_fcb(bits(28:40));
  
  % decode signs of fixed codebook first subframe, S1
  pkt.sf{1}.sPulse = decode_signs(bits(41:44));
  
  % decode gain codebook stage 1 (GA1) and stage 2 (GB1)
  pkt.sf{1}.gainIdx = decode_gain(bits(45:51));
   
  % decode second subframe 
  % decode pitch delay P2, which is relative to int(P1) - 5 2/3.  
  % 17/3 = 5 2/3
  deltaT = bin2dec(bits(52:56));
  T3 = T3 - 17 + deltaT;
  T_int = floor(T3/3);
  T_frac = mod(T3,3);
  if( T_frac == 2)
    T_int += 1;
    T_frac -= 3;
  end     

  pkt.sf{2}.T_int = T_int;
  pkt.sf{2}.T_frac = T_frac;
    
  % decode fixed codebook second subframe, C2
  pkt.sf{2}.iPulse = decode_fcb(bits(57:69));  

  % decode signs of fixed codebook second subframe, S2
  pkt.sf{2}.sPulse = decode_signs(bits(70:73));

  % decode gain codebook stage 1 (GA2) and stage 2 (GB2) 
  pkt.sf{2}.gainIdx = decode_gain(bits(74:80));
end


function sPulse = decode_signs(S)

  % convert to [s0 s1 s2 s3]
  sBits = fliplr(S);
  sPulse = arrayfun(@bin2dec,sBits);
  sPulse(sPulse == 0) = -1;    % bit = 0 decodes as sign=-1
end

function gainIdx = decode_gain(G)

  gainIdx.iA = bin2dec(G(1:3));
  gainIdx.iB = bin2dec(G(4:7));
end


function idx = decode_lsf(L)
 
  idx.i0 = bin2dec(L(1));
  idx.i1 = bin2dec(L(2:8));
  idx.i2 = bin2dec(L(9:13));
  idx.i3 = bin2dec(L(14:18));
end

function iPulse = decode_fcb(C)

  iPulse = zeros(1,4);
  i3 = bin2dec(C(1:3));
  jx = bin2dec(C(4));
  i2 = bin2dec(C(5:7));
  i1 = bin2dec(C(8:10));
  i0 = bin2dec(C(11:13));
  
  iPulse(0+(1)) = 5*i0;
  iPulse(1+(1)) = 5*i1 + 1;
  iPulse(2+(1)) = 5*i2 + 2;
  tmp = 5*i3 + 3;
  tmp += jx;
  iPulse(3+(1)) = tmp;
end
