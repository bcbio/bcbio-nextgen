function is_snp(ref, alt)
   return ref:len() == 1 and alt[1]:len() == 1
end

function is_edit_flag(ref, alt)
   if is_snp(ref, alt) then return true
   end
   return false
end
