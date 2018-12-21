package gui;

import normalizationData.AbstractNormalizationData;

public interface NormalizationDataObserver {
    public abstract void update(AbstractNormalizationData normalizationData);
}
