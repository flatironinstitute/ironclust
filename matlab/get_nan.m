function val = get_nan(S, name)
try
    val = S.(name);
catch
    val = nan;
end
end %func