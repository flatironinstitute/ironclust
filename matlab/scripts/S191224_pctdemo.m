% This function defines the outer loop that calls each generation, without
% doing the display.
function grid=callUpdateGrid(grid, gridSize, N)
    for gen = 1:N
        grid = updateGrid(grid, gridSize);
    end
end

gpuInitialGrid = gpuArray(initialGrid);

% Retain this result to verify the correctness of each version below.
expectedResult = callUpdateGrid(gpuInitialGrid, gridSize, numGenerations);

gpuBuiltinsTime = gputimeit(@() callUpdateGrid(gpuInitialGrid, ...
                                               gridSize, numGenerations));

fprintf('Average time on the GPU: %2.3fms per generation \n', ...
        1000*gpuBuiltinsTime/numGenerations);
    
    % Calculate the output value using the MEX file with shared memory. The
% initial input value is copied to the GPU inside the MEX file.
grid = pctdemo_life_mex_shmem(initialGrid, numGenerations);
gpuMexTime = gputimeit(@()pctdemo_life_mex_shmem(initialGrid, ...
                                                 numGenerations));
% Print out the average computation time and check the result is unchanged.
fprintf('Average time of %2.3fms per generation (%1.1fx faster).\n', ...
        1000*gpuMexTime/numGenerations, gpuBuiltinsTime/gpuMexTime);
assert(isequal(grid, expectedResult));


% Calculate the output value using the MEX file with textures.
grid = pctdemo_life_mex_texture(initialGrid, numGenerations);
gpuTexMexTime = gputimeit(@()pctdemo_life_mex_texture(initialGrid, ...
                                                  numGenerations));
% Print out the average computation time and check the result is unchanged.
fprintf('Average time of %2.3fms per generation (%1.1fx faster).\n', ...
        1000*gpuTexMexTime/numGenerations, gpuBuiltinsTime/gpuTexMexTime);
assert(isequal(grid, expectedResult));


fprintf('First version using gpuArray:  %2.3fms per generation.\n', ...
        1000*gpuBuiltinsTime/numGenerations);
fprintf(['MEX with shared memory: %2.3fms per generation ',...
         '(%1.1fx faster).\n'], 1000*gpuMexTime/numGenerations, ...
        gpuBuiltinsTime/gpuMexTime);
fprintf(['MEX with texture memory: %2.3fms per generation '...
         '(%1.1fx faster).\n']', 1000*gpuTexMexTime/numGenerations, ...
        gpuBuiltinsTime/gpuTexMexTime);
