function x = morse_decode(sir)
  c = morse();
  x = '*';
  for i=1:length(sir) %iau fiecare caracter din sir;
    if(sir(i) == '.') %daca este "." merg pe ramura din stanga
      c = c{2};
    elseif(sir(i) == "-")
      c = c{3};       %iar daca este "-" merg pe ramura din dreapta
    endif 
    if(length(c)==0)  %pana cand dau de ceva vid
    return;
    endif
  endfor
  x = c{1};
endfunction