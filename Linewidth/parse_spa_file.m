function S = parse_spa_file(filePath)
%PARSE_ANRITSU_SPA Parse Anritsu Spectrum Master .SPA export (text)
% Returns:
%   S.SetupInfo            (struct; grouped by sections)
%   S.TraceDataFrequency   (Nx1 double, in Hz)
%   S.TraceDataPower       (Nx1 double, in dBm)
% Also includes:
%   S.TraceFrequencyUnit   (e.g. "MHz")
%   S.TraceDataFrequencyRaw (Nx1 double, in original unit)

txt = fileread(filePath);
lines = regexp(txt, '\r\n|\n|\r', 'split');
n = numel(lines);

% -------------------- Setup parsing (everything before "# Begin SPA Data") --------------------
iBeginData = find(contains(lower(strtrim(lines)), 'begin spa data'), 1, 'first');
if isempty(iBeginData); iBeginData = n+1; end

SetupInfo = struct();
section = "Header";  % default bucket

for i = 1:min(iBeginData-1, n)
    L = strtrim(lines{i});
    if L == "" || startsWith(L, "<") || startsWith(L, "#") || startsWith(L, ">")
        % Capture section changes like "# Begin TRACE A Setup"
        if startsWith(L, "# Begin") && contains(L, "Setup")
            name = regexprep(L, '^#\s*Begin\s*', '');
            name = regexprep(name, '\s*Setup\s*$', '');
            name = strrep(strtrim(name), " ", "");
            section = string(matlab.lang.makeValidName(name));
            if ~isfield(SetupInfo, section), SetupInfo.(section) = struct(); end
        end
        continue
    end

    % Expect "KEY,VALUE"
    parts = strsplit(L, ',');
    if numel(parts) < 2, continue; end
    key = matlab.lang.makeValidName(strtrim(parts{1}));
    valStr = strtrim(strjoin(parts(2:end), ',')); % in case value contains commas

    num = str2double(valStr);
    if ~isnan(num) && ~isempty(valStr)
        val = num;
    else
        val = valStr;
    end

    if ~isfield(SetupInfo, section), SetupInfo.(section) = struct(); end
    SetupInfo.(section).(key) = val;
end

% -------------------- Trace parsing (P_# lines) --------------------
pMask = false(n,1);
for i = 1:n
    pMask(i) = ~isempty(regexp(strtrim(lines{i}), '^P_\d+\s*,', 'once'));
end
pLines = lines(pMask);

N = numel(pLines);
pow = nan(N,1);
freqRaw = nan(N,1);
unit = "";

for k = 1:N
    parts = strtrim(strsplit(pLines{k}, ','));
    % Typical: P_0, <power> , <freq> ,MHz
    if numel(parts) >= 3
        pow(k) = str2double(parts{2});
        freqRaw(k) = str2double(parts{3});
    end
    if numel(parts) >= 4 && unit == ""
        unit = string(parts{4});
    end
end

mult = unitMultiplier(unit);
freqHz = freqRaw * mult;

S = struct();
S.SetupInfo = SetupInfo;
S.TraceDataFrequency = freqHz;
S.TraceDataPower = pow;
S.TraceFrequencyUnit = unit;
S.TraceDataFrequencyRaw = freqRaw;
end

function mult = unitMultiplier(u)
u = upper(string(strtrim(u)));
switch u
    case "HZ",  mult = 1;
    case "KHZ", mult = 1e3;
    case "MHZ", mult = 1e6;
    case "GHZ", mult = 1e9;
    otherwise,  mult = 1; % fallback
end
end
