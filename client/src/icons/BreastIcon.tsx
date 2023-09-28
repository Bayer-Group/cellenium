import Icon from './svg/breast.svg';

function BreastIcon({ size }: { size: number }) {
  return (
    <div style={{ width: size }}>
      <img src={Icon} alt="breast icon" />
    </div>
  );
}

export default BreastIcon;
