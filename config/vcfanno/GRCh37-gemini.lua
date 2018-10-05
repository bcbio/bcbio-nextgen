function mean(vals)
    local sum=0
	for i=1,#vals do
		sum = sum + vals[i]
	end
	return sum / #vals
end

function loc(chrom, start, stop)
    return chrom .. ":" .. start .. "-" .. stop
end

CLINVAR_LOOKUP = {}
CLINVAR_LOOKUP['0'] = 'unknown'
CLINVAR_LOOKUP['1'] = 'germline'
CLINVAR_LOOKUP['2'] = 'somatic'
CLINVAR_LOOKUP['4'] = 'inherited'
CLINVAR_LOOKUP['8'] = 'paternal'
CLINVAR_LOOKUP['16'] = 'maternal'
CLINVAR_LOOKUP['32'] = 'de-novo'
CLINVAR_LOOKUP['64'] = 'biparental'
CLINVAR_LOOKUP['128'] = 'uniparental'
CLINVAR_LOOKUP['256'] = 'not-tested'
CLINVAR_LOOKUP['512'] = 'tested-inconclusive'
CLINVAR_LOOKUP['1073741824'] = 'other'

CLINVAR_SIG = {}
CLINVAR_SIG['0'] = 'uncertain'
CLINVAR_SIG['1'] = 'not-provided'
CLINVAR_SIG['2'] = 'benign'
CLINVAR_SIG['3'] = 'likely-benign'
CLINVAR_SIG['4'] = 'likely-pathogenic'
CLINVAR_SIG['5'] = 'pathogenic'
CLINVAR_SIG['6'] = 'drug-response'
CLINVAR_SIG['7'] = 'histocompatibility'
CLINVAR_SIG['255'] = 'other'
CLINVAR_SIG['.'] = '.'

function intotbl(ud)
    local tbl = {}
    for i=1,#ud do
        tbl[i] = ud[i]
    end
    return tbl
end

-- from lua-users wiki
function split(str, sep)
        local sep, fields = sep or ":", {}
        local pattern = string.format("([^%s]+)", sep)
        str:gsub(pattern, function(c) fields[#fields+1] = c end)
        return fields
end

function contains(str, tok)
    return string.find(str, tok) ~= nil
end


function div2(a, b)
    if(a == 0) then return "0.0" end
    return string.format("%.9f", (a + 0) / b)
end

function ratio(vals)
    vals = vals[1] -- get 2 values per element. ref and alt counts.
    if vals[2] == 0 then return "0.0" end
    return string.format("%.9f", vals[2] / (vals[1] + vals[2]))
end

function clinvar_sig(vals)
    local t = type(vals)
    -- just a single-value
    if(t == "string" or t == "number") and not contains(vals, "|") then
        return CLINVAR_SIG[vals]
    elseif t ~= "table" then
        if not contains(t, "userdata") then
            if t == "string" then
                vals = split(vals, ",")
            else
                vals = {vals}
            end
        else
            vals = intotbl(vals)
        end
    end
    local ret = {}
    for i=1,#vals do
        if not contains(vals[i], "|") then
            ret[#ret+1] = CLINVAR_SIG[vals[i]]
        else
            local invals = split(vals[i], "|")
            local inret = {}
            for j=1,#invals do
                inret[#inret+1] = CLINVAR_SIG[invals[j]]
            end
            ret[#ret+1] = join(inret, "|")
        end
    end
    return join(ret, ",")
end

join = table.concat

function check_clinvar_aaf(clinvar_sig, max_aaf_all, aaf_cutoff)
    -- didn't find an aaf for this so can't be common
    if max_aaf_all == nil or clinvar_sig == nil then
        return false
    end
    if type(clinvar_sig) ~= "string" then
        clinvar_sig = join(clinvar_sig, ",")
    end
    if false == contains(clinvar_sig, "pathogenic") then
        return false
    end
    if type(max_aaf_all) ~= "table" then
        return max_aaf_all > aaf_cutoff
    end
    for i, aaf in pairs(max_aaf_all) do
        if aaf > aaf_cutoff then
            return true
        end
    end
    return false
end

function check_population_aaf(max_aaf_all, aaf_cutoff)
    -- didn't find an aaf for this so can't be common
    if max_aaf_all == nil then
        return false
    end
    if type(max_aaf_all) ~= "table" then
        return max_aaf_all > aaf_cutoff
    end
    for i, aaf in pairs(max_aaf_all) do
        if aaf > aaf_cutoff then
            return true
        end
    end
    return false
end


function setid(...)
	local t = {...}
	local res = {}
	local seen = {}
	for i, v in pairs(t) do
		if v ~= "." and v ~= nil and v ~= "" then
			if seen[v] == nil then
				res[#res+1] = string.gsub(v, ",", ";")
				seen[v] = true
			end
		end
	end
	return table.concat(res, ";")
end
