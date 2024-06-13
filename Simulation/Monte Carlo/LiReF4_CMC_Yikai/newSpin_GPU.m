function spins = newSpin_GPU(ion, basis, wav, spinType)
    % Get the correct spin matrices based on type
    switch spinType
        case 'electron'
            Jx = gpuArray(ion.Jx);
            Jy = gpuArray(ion.Jy);
            Jz = gpuArray(ion.Jz);
        case 'nuclear'
            Jx = gpuArray(ion.Ix);
            Jy = gpuArray(ion.Iy);
            Jz = gpuArray(ion.Iz);
        otherwise
            error('Unrecognized spin type!');
    end
    % Compute spins for all wavefunctions
    spx = real(wav' * basis' * Jx * basis * wav);
    spy = real(wav' * basis' * Jy * basis * wav);
    spz = real(wav' * basis' * Jz * basis * wav);
    spins = gpuArray([spx; spy; spz]');
end
