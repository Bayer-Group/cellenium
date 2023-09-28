import Icon from './svg/esophagus.svg';

function EsophagusIcon({ size }: { size: number }) {
  return (
    <div style={{ width: size }}>
      <img src={Icon} alt="esophagus icon" />
    </div>
  );
}

export default EsophagusIcon;
