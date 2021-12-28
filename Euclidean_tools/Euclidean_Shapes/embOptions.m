function opts = embOptions(opts_custom)

% set default options
opts = struct('ThetaSpline','complete',...
              'AffineTrans','LA',...
              'TE','sharp',...
              'M',[],...
              'b',[],...
              'Minv',[]);
          
if ~isempty(opts_custom)
    all_fields = fieldnames(opts);
    for field = fieldnames(opts_custom)'
        tf_fieldmatch = strcmpi(all_fields, field);
        if any(tf_fieldmatch)
            opts = setfield(opts, all_fields{tf_fieldmatch}, getfield(opts_custom, field{:}));
        else
            error(['ERROR: Unexpected field name -> ' field])
        end
    end
end

end