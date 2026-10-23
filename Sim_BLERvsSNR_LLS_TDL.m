function Sim_BLERvsSNR_LLS_TDL(SNR_list,...
    MCS,TBS,CR,modulation,...
    NSizeGrid,Nfft, channel_type,fd, DS,...
    N, W, sigma_2)
% Link-level Simulation for FAS
% HJ Hong
% BLER vs SNR 
% Input:
        % SNR_list
% MCS:
% TBS:
% CR: Target Code Rate x 1024
% modulation: 
% NSizeGrid: number of grid
% channel_type: 'TDL-A'/'TDL-B'/'TDL-C'/'TDL-D'/'TDL-E' 
% fd: Doppler frequency
% DS: Delay Spread
% N: Number of antenna ports in FAS
% W: normalized physical size of FAS
% sigma_2: $\sigma^2$ the power of the channel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    NN_block_max = 20000;
    NN_block_par = 100;
    dir_name = ['result_tmp_5MHz/',datestr(now,30),'_LLS_TDL_MCS',num2str(MCS),...
        '_D',num2str(fd),'Hz_DS',num2str(DS*1e9),'us'];
    % 5G NR
    Target_Code_rate = CR/1024;
    trBlkSizes = TBS;
    nlayers = 1;
    rv = 0;
    N_in = 50;
    % fluid antenna
    N_port = N(1)*N(2);% number of fluid antenna ports  
    FAMA_channel_type = 4;
    %%%%%%%%%%%%%%%% Initiating %%%%%%%%%%%%%%%%%%%%%%%%%%%
    BER_list = zeros(size(SNR_list));
    BLER_list = zeros(size(SNR_list));
    %%%%%%%%%%%%% iterative simulation %%%%%%%%%%%%%%%%%%%%
    for idx = 1:length(SNR_list)
        %%%%%%%%% initialization %%%%%%%%%%%%%%%%%%%%
        SNR = SNR_list(idx);
        snr = 10^(SNR/10);
        N0 = sigma_2/(snr*double(Nfft));
        N_block = 0;
        Num_Blk_error = 0;
        Num_Bit_error = 0;
        %%%%%%%% N_block simulations %%%%%%%%%%%%%%%%
        while (N_block < NN_block_max)
            %%%% initialization %%%%%%%%%%%%%%
            Blk_error_par = zeros(1,NN_block_par);
            Bit_error_par = zeros(1,NN_block_par);
            parfor idx_block = 1:NN_block_par
                %%%% channel %%%%%%%%%%%%%%%
                [channel, maxChDelay] = getTDLchannel(NSizeGrid, channel_type, ...
                    N, W, sigma_2, FAMA_channel_type, fd, DS);
                [Blk_error_par(idx_block), Bit_error_par(idx_block)] = Sim1Frame(...
                    trBlkSizes, modulation, Target_Code_rate, ...
                    NSizeGrid, nlayers, rv, N_in, channel, maxChDelay, N0,1/snr);
                if mod(idx_block,50)==0
                    fprintf('SNR: %.2f, Block test: %d.\n',...
                        SNR, N_block + idx_block)
                end
            end
            N_block = N_block +NN_block_par;
            Num_Blk_error = Num_Blk_error +sum(Blk_error_par);
            Num_Bit_error = Num_Bit_error +sum(Bit_error_par);
            fprintf(['SNR: %.2f, Block test: %d, Block error: %d, Bit error: %d,\n'],...
                SNR, N_block, Num_Blk_error, Num_Bit_error)
            if ~exist(dir_name,'dir')
                mkdir(dir_name)
            end
            save([dir_name,'/result_tmp_MCS',num2str(MCS),'_SNR',num2str(SNR),'.mat'],...
                'N_block','SNR','trBlkSizes',...
                'Num_Bit_error','Num_Blk_error')
            %%% conditions of break while:
            if (Num_Blk_error >= 20) && (mod(N_block,50) == 0)
                break
            elseif (Num_Blk_error >= 10) &&(mod(N_block,100) == 0)
                break
            elseif (Num_Blk_error >= 5) &&(mod(N_block,500) == 0)
                break
            elseif (Num_Blk_error >= 2) &&(mod(N_block,1000) == 0)
                break
            elseif (Num_Blk_error >= 1) &&(mod(N_block,5000) == 0)
                break
            end
        end
        BER_list(idx) = Num_Bit_error/TBS/N_block;
        BLER_list(idx) = Num_Blk_error/N_block;
        save([dir_name,'/SNRvsBLER_MCS',num2str(MCS),'_W',num2str(W(1)),'x',num2str(W(2)),'_N',num2str(N(1)),'x',num2str(N(2)),'.mat'],...
            'SNR_list','BLER_list','BER_list')
        if (Num_Blk_error/N_block<=5e-3)
            break
        end
    end
    fprintf(['\n',datestr(now,31),'\n'])
end

function [channel, maxChDelay] = getTDLchannel(NSizeGrid, channel_type, ...
    N, W, sigma_2, FAMA_channel_type, fd, DS)
% channel_type: "TDL-A", "TDL-B", "TDL-C", "TDL-D" or "TDL-E"
    N_port = N(1)*N(2); 
    G = channel_g_FAS(N,W,sigma_2,FAMA_channel_type);
    carrier = nrCarrierConfig;
    carrier.NSizeGrid = NSizeGrid;
    ofdmInfo = nrOFDMInfo(carrier);
    channel = nrTDLChannel;
    channel.DelayProfile = channel_type; 
    channel.SampleRate = ofdmInfo.SampleRate;
    channel.NumTransmitAntennas = 1;
    channel.NumReceiveAntennas = N_port;
    channel.MIMOCorrelation = 'Custom';
    if N_port > 1
        channel.ReceiveCorrelationMatrix = G.SSS;
    else
        channel.ReceiveCorrelationMatrix = 1;
    end
    channel.NormalizeChannelOutputs = 0;
    channel.TransmitCorrelationMatrix = 1;
    channel.Seed = randi(intmax,1);
    channel.MaximumDopplerShift = fd;
    channel.DelaySpread = DS;
    chInfo = info(channel);
    maxChDelay = ceil(max(chInfo.PathDelays*channel.SampleRate)) + chInfo.ChannelFilterDelay;
end

function [Blk_error, Bit_error] = Sim1Frame(...
    trBlkSizes, modulation, TargetCoderate, ...
    NSizeGrid, nlayers, rv, N_in, ...
    channel,maxChDelay, N0_TD,N0_FD)
% BLER Simulation for the 1 subframe
% Input:
%   trBlkSizes: TBS
%   modulation: 'QPSK', '16QAM', '64QAM'
%   TargetCodeRate: Target Code Rate
%   NSizeGrid: N_RB, number of PRB
%   nlayers: 1
%   rv: 0
%   N_in: number of decoding iterations
% Output:
%   Blk_error_List: Block error rate
%   Bit_error_List: bit error rate
    %%%% transmit %%%%%%%%%%%%%%%%%
    data = randi([0 1],trBlkSizes,1,'int8');
    [TransSig, pdschSym,~] = Transmitter_LLS(data,...
        NSizeGrid, modulation, nlayers, TargetCoderate);
    TransSig = [TransSig; zeros(maxChDelay,size(TransSig,2))];
    [rxWaveform,pathGains,sampleTimes] = channel(TransSig);
    nn = sqrt(N0_TD/2)*(randn(size(rxWaveform))+1j*randn(size(rxWaveform)));
    r = rxWaveform+nn;
    %%%% Receiver %%%%%%%%%%%%%%%%%%%%
    % time sync.
    pathFilters = getPathFilters(channel);
    [offset, ~] = nrPerfectTimingEstimate(pathGains, pathFilters);
    r = r(1+offset:end,:);
    % r_tmp = r(:,PortSeq);
    % OFDM Dem.
    [rxGrid] = OFDMDem(r.', NSizeGrid, size(r,2));
    % channel Est.
    estChGridAnts = PerfectChanEst(pathGains, pathFilters, offset, sampleTimes, NSizeGrid);
    % Equalization
    r_pdsch = ExtractSym(rxGrid, NSizeGrid, modulation, nlayers);
    g_perfect = ExtractSym(estChGridAnts, NSizeGrid, modulation, nlayers);
    SNR_k = mean(abs(g_perfect).^2,2)./N0_FD;
    [~, k_star] = max(SNR_k);
    r_s = r_pdsch(k_star,:);
    g_s = g_perfect(k_star,:);
    W_MMSE = conj(g_s)./(conj(g_s).*g_s+N0_FD);
    r_MMSE = W_MMSE.*r_s;
    [dlschLLRs,~] = Dem(r_MMSE, N0_FD*mean(abs(W_MMSE).^2), NSizeGrid, modulation, nlayers);
    [data_out_AP, ~] = DLSCHDec(dlschLLRs, TargetCoderate, trBlkSizes, rv, N_in, modulation, nlayers);
    Blk_error = (sum(data~=data_out_AP) > 0);
    Bit_error = (sum(data~=data_out_AP));
end

function [TransSig, pdschSym, waveformInfo] = Transmitter_LLS(InputBits,...
    NSizeGrid, modulation, nlayers, Target_Code_rate)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data Process at the transmitter
% LLS
% ----------------------------------------------
% Input:
% InputBits: Bits to be transmitted, TBS x U bits
% encodeDLSCH: DLSCH encoder
% NSizeGrid: ABW/180kHz
% ---------------------------------------------
% Output:
% TransSig: Transmitted signal,  x U symbols
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    carrier = nrCarrierConfig;
    carrier.NSizeGrid = NSizeGrid;
    pdsch = nrPDSCHConfig;
    pdsch.PRBSet = 0:NSizeGrid-1;
    pdsch.Modulation = modulation;
    pdsch.NumLayers = nlayers;
    encodeDLSCH = nrDLSCH;
    encodeDLSCH.MultipleHARQProcesses = false;
    encodeDLSCH.TargetCodeRate = Target_Code_rate;
    [pdschIndices,pdschInfo] = nrPDSCHIndices(carrier,pdsch);
    trBlk = InputBits;
    encodeDLSCH.release
    setTransportBlock(encodeDLSCH,trBlk,0,0);
    codedTrBlock = encodeDLSCH(pdsch.Modulation,pdsch.NumLayers,pdschInfo.G,0);
        % carrier.NCellID = idx_UU;
    pdschSym = nrPDSCH(carrier,pdsch,codedTrBlock);
        % DMRS
    dmrsSymbols = nrPDSCHDMRS(carrier,pdsch);
    dmrsIndices = nrPDSCHDMRSIndices(carrier,pdsch);
    pdschGrid = nrResourceGrid(carrier,1);
    [~,pdschAntIndices] = nrExtractResources(pdschIndices,pdschGrid);
    pdschGrid(pdschAntIndices) = pdschSym;
    [~,dmrsAntIndices] = nrExtractResources(dmrsIndices,pdschGrid);
    pdschGrid(dmrsAntIndices) = pdschGrid(dmrsAntIndices) + dmrsSymbols;
    [TransSig,waveformInfo] = nrOFDMModulate(carrier,pdschGrid);
end
function [rxGrid] = OFDMDem(RecSig, NSizeGrid, N_port)
    carrier = nrCarrierConfig;
    carrier.NSizeGrid = NSizeGrid;
    rxGrid = zeros(carrier.NSizeGrid*12,14,N_port);
    for idx_port = 1:N_port
        rxGrid(:,:,idx_port) = nrOFDMDemodulate(carrier,RecSig(idx_port,:).');
    end
end
function [ExtSym] = ExtractSym(rxGrid, NSizeGrid, modulation, nlayers)
    carrier = nrCarrierConfig;
    carrier.NSizeGrid = NSizeGrid;
    pdsch = nrPDSCHConfig;
    pdsch.PRBSet = 0:NSizeGrid-1;
    pdsch.Modulation = modulation;
    pdsch.NumLayers = nlayers;
    [pdschIndices,~] = nrPDSCHIndices(carrier,pdsch);
    ExtSym = nrExtractResources(pdschIndices,rxGrid).';
end
function [estChGridAnts] = PerfectChanEst(pathGains, pathFilters, offset, ...
    sampleTimes, NSizeGrid)
    carrier = nrCarrierConfig;
    carrier.NSizeGrid = NSizeGrid;
    estChGridAnts = nrPerfectChannelEstimate(carrier,pathGains,pathFilters,offset,sampleTimes);
end

function [dlschLLRs, rxsyms] = Dem(r_MMSE, N0, ...
    NSizeGrid, modulation, nlayers)
    carrier = nrCarrierConfig;
    carrier.NSizeGrid = NSizeGrid;
    pdsch = nrPDSCHConfig;
    pdsch.PRBSet = 0:NSizeGrid-1;
    pdsch.Modulation = modulation;
    pdsch.NumLayers = nlayers;
    [dlschLLRs, rxsyms] = nrPDSCHDecode(carrier,pdsch,r_MMSE.',N0);
end

function [data_out, blkerr] = DLSCHDec(dlschLLRs, ...
    Target_Code_rate, trBlkSizes, rv, N_in, modulation, nlayers)
    decodeDLSCH = nrDLSCHDecoder;
    decodeDLSCH.MultipleHARQProcesses = false;
    decodeDLSCH.TargetCodeRate = Target_Code_rate;
    decodeDLSCH.LDPCDecodingAlgorithm = "Normalized min-sum";
    decodeDLSCH.MaximumLDPCIterationCount = N_in;
    decodeDLSCH.TransportBlockLength = trBlkSizes;
    [data_out,blkerr] = decodeDLSCH(dlschLLRs,modulation,nlayers,rv);
end