function prot_info = getProtInfoUniProt(acc, from, to, fields)
%% prot_name = getProtInfoUniProt(acc, from, to, columns)
% Use the UniProt REST API to retrieve information associated with UniProt
% entries for a list of protein IDs.
%
% Input
%   cellstr acc:       list of IDs that should be mapped or for which
%                      information should be retrieved
%   char from:         name of database, where IDs in acc input originate
%                      from (by default 'UniProtKB_AC-ID')
%   char to:           name of database, to which IDs should be translated
%                      (by default: 'UniProtKB')
%   cellstr fields:   array of column names that should be returned
% Output
%   table prot_info:    table containing the requested protein information

% warning('This function currently returns the whole result structure as result entry.')
% The reason is that in response to given fields, the API returns multiple
% results.

if ~ischar(acc) && ~iscellstr(acc)
    cell_idx = ~cellfun(@ischar, acc);
    acc(cell_idx) = acc{cell_idx};
elseif ischar(acc)
    acc = cellstr(acc);
end

if nargin >= 2 && ~isempty(from)
    if ~ischar(from)
        error('""from"" input must be of type char')
    end
else
    from = 'UniProtKB_AC-ID';
end

if nargin >= 3 && ~isempty(to)
    if ~ischar(to)
        error('""to"" input must be of type char')
    end
else
    to = 'UniProtKB';
end

if nargin == 4
    if ischar(fields)
        fields = cellstr(fields);
    elseif ~iscellstr(fields)
        cell_idx = ~cellfun(@ischar, fields);
        fields(cell_idx) = fields{cell_idx};
    end
else
    fields = {'accession'};
end

URL = 'https://rest.uniprot.org/idmapping/run';

% divide gene names into chunks
chunk_size = 100;
startIdx = 1:chunk_size:numel(acc);
prot_info = table;
for i=1:numel(startIdx)
    if i<numel(startIdx)
        endIdx = startIdx(i)+(chunk_size-1);
    else
        endIdx = numel(acc);
    end
    fprintf('Processing IDs %d to %d ...\n', startIdx(i), endIdx)
    tmpIdx = startIdx(i):endIdx;
    tmpQuery = acc(tmpIdx);
    
    curl_cmd_jobid = ['curl --silent --request POST ' URL ' ',...
        '--form ids="', strjoin(tmpQuery,',') '" ',...
        '--form from="' from '" ',...
        '--form to="' to '"'];
    % send requestto retrieve job ID
    status = 1;
    trials = 0;
    while status ~= 0 && trials < 3
        trials = trials + 1;
        [status, job_id_json] = system(curl_cmd_jobid);
    end
    
    if status ~= 0 || ~isfield(jsondecode(job_id_json), 'jobId')
        error('Request not successful after %i attempts', trials)
    else
        job_id = jsondecode(job_id_json).jobId;
    end
    
    % wait until job is finished
    job_status = '';
    curl_cmd_jobstatus = ['curl --silent https://rest.uniprot.org/idmapping/status/' ...
        job_id];
    tic
    t = 0;
    while ~isequal(job_status, 'FINISHED') && t < 120
        t = toc;
        [status, job_status_json] = system(curl_cmd_jobstatus);
        if status == 0
            job_status = jsondecode(job_status_json).jobStatus;
        end
    end
    
    % get job results
    curl_cmd_jobresult = ['curl --silent https://rest.uniprot.org/idmapping/uniprotkb/results/stream/' ...
        job_id '?fields=' strjoin(fields, ',')];
    [status, job_result_json] = system(curl_cmd_jobresult);
    if status == 0
        job_result_struct = jsondecode(job_result_json).results;
        job_result_table = struct2table(job_result_struct);
    else
        error('Job results could not be fetched')
    end
    prot_info = [prot_info; job_result_table];
end

end