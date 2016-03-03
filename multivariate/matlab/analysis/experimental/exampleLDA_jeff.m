clear all;
tic;

%% Load corpus
corpus = 'wiki'; %'wiki';%'tasa'
modifier = '';%'_large';%'_reduced_mturk';%
load(sprintf('%s/corpus_%s%s.mat', corpus,corpus, modifier));
doc_indices = double(doc_indices);

%% Get a stopword list
stopwords = textread( 'stopwordlist.txt' , '%s' );

%%
% Set the number of topics
T=1000;%300;%200; % Changed by DC

%%
% What output to show (0=no output; 1=iterations; 2=all output)
OUTPUT = 1;

%%
% The number of iterations
BURNIN   = 1000; % the number of iterations before taking samples (change to 1000 for real runs!)

%% Filter out the stopwords
% [ ismem , badlist ] = ismember( stopwords , vocabulary );
% badlist = badlist( ismem == 1 );
% okplaces = find( ~ismember( word_indices , badlist ));
% word_indices = word_indices( okplaces );
% doc_indices = doc_indices( okplaces );
% doc_indices = double(doc_indices); % Added by DC

%% Moved and changed by DC
% Set the hyperparameters
BETA = 200 / length(vocabulary);%0.01;
ALPHA = 50 / T;%0.1;

%%
% This is to run just 1 chain and collect 1 sample.
% N = BURNIN;
% fprintf( 'Running Gibbs sampler for burnin\n' );
% [ WP,DP,Z ] = GibbsSamplerLDA( wordindex , docindex , T , N , ALPHA , BETA , SEED , OUTPUT );
% wp = WP;
% save(sprintf('tasa300/wp%d.mat', SEED), 'wp');

% This is to run multiple chains, collecting multiple samples from each
% chain. The outputs from Mark Steyvers seem to be generated this way.
LAG      = 100; % the lag between samples
NSAMPLES = 8; % the number of samples for each chain
NCHAINS  = 1;%3;% % the number of chains to run

SEED = 0;%8;%16;

for chain = 1:NCHAINS
    fprintf( sprintf('Running Gibbs sampler for burnin, chain %d, sample 1\n', chain) );
    tic;
    SEED = SEED + 1;
    N = BURNIN;
    [ wp,DP,Z ] = GibbsSamplerLDA( word_indices , doc_indices , T , N , ALPHA , BETA , SEED , OUTPUT );
    save(sprintf('%s/%d/wp%s_chain%d_sample1.mat', corpus, T, modifier, chain), 'wp');   % Added by DC
    toc;
    
    fprintf( 'Continue to run sampler to collevvct samples\n' );
    for sample = 2:NSAMPLES
        fprintf('Sample %d\n', sample);
        tic;
        N = LAG;
        SEED = SEED + 1; % important -- change the seed between samples !!
        [ wp,DP,Z ] = GibbsSamplerLDA( word_indices , doc_indices , T , N , ALPHA , BETA , SEED , OUTPUT , Z );
        save(sprintf('%s/%d/wp%s_chain%d_sample%d.mat', corpus, T, modifier, chain, sample), 'wp');  % Added by DC
        toc;
    end
end

%topicWords = WriteTopics( WP , BETA , vocabulary , 30 , 1.0 , 4 , 'tasa300/topic_words_top30.txt' ); % Changed by DC

%toc;



