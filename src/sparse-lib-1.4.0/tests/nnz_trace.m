function nnz_trace( fn, m, n, step );
    %Compile sbd2trp with the _SBD2TRP_COUT flag,
    %and without the _INFO or _DEBUG flags.
    %Run sbd2trp on any m times n matrix,
    %and pipe the stdout to a file (e.g. 'trace').
    %Then run "nnz_trace('trace',<m>,<n>)" (with
    %the appropiate values for <m> and <n>) to
    %visually inspect the induced ordering.
    %Pressing any key in the main MATLAB window
    %will proceed to the next nonzero.
    %Hold a key (e.g., space bar) for a 'movie'
    %effect.
    if( nargin < 4 )
        step = 1;
    end

    A=dlmread(fn);
    figure(1);clf;
    plot(n,m); hold on;
    plot(0,0); hold on;
    c=0;
    for a=A'
        plot(a(2),m-a(1),'.'); hold on;
        c=c+1;
        if mod(c,step)==0
            pause;
        end
    end
end
